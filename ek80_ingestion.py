"""EK80 ingestion and calibration.

Reads raw Simrad EK80 ``.raw`` files with Echopype, writes a lossless
SONAR-netCDF4 archive copy, and calibrates the complex backscatter to
volume backscattering strength (Sv, dB re 1 m-1).

Key interfaces:
    load_config()      read cruises/<cruise_id>/config.toml
    ingest_ek80_raw()  .raw -> EchoData (+ optional .nc archive)
    calibrate_sv()     EchoData -> calibrated Sv Dataset
"""

import logging
import pathlib
import re
import struct
import tomllib

import echopype as ep
import numpy as np
import pandas as pd
import xarray as xr
from echopype.utils.uwa import calc_sound_speed

REQUIRED_SECTIONS = ('environment', 'platform', 'processing', 'output')


def _flatten_toml(raw, path):
    flat = {}
    for section, values in raw.items():
        if not isinstance(values, dict):
            flat[section] = values
            continue
        for key, value in values.items():
            if key in flat:
                raise ValueError(f"Duplicate config key '{key}' in {path}")
            flat[key] = value
    return flat


def load_config_with_defaults(cruise_config_path, defaults_path):
    """Cruise ``config.toml`` layered over the project's ``config.defaults.toml``.

    A key set in the cruise file always wins; a missing key falls back to the defaults file.  The
    cruise file may be absent altogether.  Returns the flat config dict, which also carries
    ``_CONFIG_SOURCE`` (what was used) so that every build records it.
    """
    cruise_config_path, defaults_path = pathlib.Path(cruise_config_path), pathlib.Path(defaults_path)
    have_cruise, have_defaults = cruise_config_path.exists(), defaults_path.exists()
    if not have_cruise and not have_defaults:
        raise FileNotFoundError(
            f"No config: {cruise_config_path} and {defaults_path} are both missing. "
            f"Run 'python main.py <cruise_id> --init-config' to write a starter config.toml.")
    config = load_config(defaults_path) if have_defaults else {}
    if have_cruise:
        config.update(load_config(cruise_config_path))
        config['_CONFIG_SOURCE'] = ('cruise config.toml over config.defaults.toml' if have_defaults
                                    else 'cruise config.toml')
    else:
        config['_CONFIG_SOURCE'] = 'config.defaults.toml ONLY (the cruise has no config.toml)'
    return config


def resolve_channel_settings(config, ds_Sv):
    """Fills in the channel choices a config may leave out, from the frequencies in this file.

    The seabed is detected on, and schools are searched in, the lowest-frequency channel unless the
    config names one: low frequencies reach deepest and give the cleanest seabed echo.
    Returns a copy; the original is not changed.
    """
    config = dict(config)
    lowest = float(np.min(ds_Sv['frequency_khz'].values))
    if 'SEABED_CHANNEL_KHZ' not in config:
        config['SEABED_CHANNEL_KHZ'] = lowest
        logging.info(f"SEABED_CHANNEL_KHZ not set: using the lowest frequency in this file ({lowest:g} kHz)")
    if 'SCHOOL_FREQUENCIES_KHZ' not in config:
        config['SCHOOL_FREQUENCIES_KHZ'] = [lowest]
    return config


def load_config(path):
    """Loads a cruise ``config.toml`` and returns it as a flat dict.

    Section names are only for readability; every key must be unique across
    sections.  Optional keys that are absent are simply not in the dict.
    """
    path = pathlib.Path(path)
    with open(path, 'rb') as f:
        raw = tomllib.load(f)
    return _flatten_toml(raw, path)


def _nmea_degrees(value, hemisphere):
    """NMEA ddmm.mmmm / dddmm.mmmm plus N/S/E/W -> signed decimal degrees."""
    deg = int(float(value) // 100)
    out = deg + (float(value) - deg * 100) / 60.0
    return -out if hemisphere in ('S', 'W') else out


def read_gps_fixes(path):
    """GPS fixes recorded inside a Simrad .raw file (NMEA ``GGA`` sentences in its navigation datagrams).

    Returns a DataFrame with ``time`` (UTC, from the datagram's own PC timestamp: the same clock as the
    pings, which can differ from the GPS clock), ``lat``, ``lon``, ``quality`` (GGA fix quality; 0 =
    invalid fixes are dropped) and ``satellites``.  Empty when the file has no usable GPS.
    Only the small navigation datagrams are read; the bulky echo datagrams are skipped.
    """
    rows = []
    with open(path, 'rb') as f:
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            length = struct.unpack('<i', head)[0]
            kind = f.read(4)
            if kind == b'NME0':
                body = kind + f.read(length - 4)
                f.read(4)
                try:
                    stamp = struct.unpack('<Q', body[4:12])[0]            # Windows FILETIME: 100 ns since 1601
                    sentence = body[12:].split(b'\x00')[0].decode('ascii', 'ignore')
                    if sentence[:1] != '$' or sentence[3:6] != 'GGA':
                        continue
                    p = sentence.split('*')[0].split(',')
                    quality = int(p[6] or 0)
                    if quality == 0 or not p[2] or not p[4]:
                        continue
                    t = pd.Timestamp(1970, 1, 1) + pd.Timedelta(microseconds=(stamp - 116444736000000000) // 10)
                    rows.append((t, _nmea_degrees(p[2], p[3]), _nmea_degrees(p[4], p[5]), quality,
                                 int(p[7]) if p[7] else np.nan))
                except (ValueError, IndexError, struct.error):
                    continue
            else:
                f.seek(length - 4 + 4, 1)                                  # skip the rest of this datagram and its length
    return pd.DataFrame(rows, columns=['time', 'lat', 'lon', 'quality', 'satellites'])


def _fmt_deg_min(value, pos, neg, width):
    hemi = pos if value >= 0 else neg
    value = abs(value)
    deg = int(value)
    return f"{deg}-{(value - deg) * 60:05.2f}{hemi}"


def _spread_m(fixes):
    """Largest north-south or east-west extent of a set of fixes, in metres."""
    ns = (fixes['lat'].max() - fixes['lat'].min()) * 111_000.0
    ew = (fixes['lon'].max() - fixes['lon'].min()) * 111_000.0 * np.cos(np.radians(fixes['lat'].mean()))
    return float(max(ns, ew))


def stations_from_gps(raw_files, moving_m=500.0):
    """One station row per survey from the GPS inside the raw files.

    ``raw_files`` is a list of paths.  Returns (rows, fixes_by_survey): rows are dicts with the
    stations.csv columns (survey_id, lat, lon, lat_decimal, lon_decimal, confidence, note); a survey
    with no usable GPS gets no row.  Surveys whose fixes spread more than ``moving_m`` are flagged in
    the note (one position is then a poor description) and the caller may write a nav.csv.
    """
    by_survey = {}
    for p in raw_files:
        fixes = read_gps_fixes(p)
        if len(fixes):
            by_survey.setdefault(survey_id_from_file_key(raw_file_key(pathlib.Path(p))), []).append(fixes)
    rows, out = [], {}
    for sid, parts in by_survey.items():
        fixes = pd.concat(parts, ignore_index=True)
        out[sid] = fixes
        la, lo = float(fixes['lat'].median()), float(fixes['lon'].median())
        spread = _spread_m(fixes)
        moving = spread > moving_m
        q = fixes['quality'].mode().iloc[0]
        confidence = 'fair' if (moving or q == 1) else 'good'
        note = (f"Written automatically from the median of {len(fixes)} GPS (NMEA GGA) fixes inside {len(parts)} raw file(s); "
                f"fix quality {int(q)} (1 = GPS, 2 = DGPS, 4 = RTK fixed, 5 = RTK float); the fixes spread {spread:.0f} m"
                + (": the ship MOVED, so one position is a poor description - use the track in nav.csv. " if moving else ". ")
                + "The GPS clock can differ from the file-name clock; ping times are matched using the file clock.")
        rows.append(dict(survey_id=sid, lat=_fmt_deg_min(la, 'N', 'S', 2), lon=_fmt_deg_min(lo, 'E', 'W', 3),
                         lat_decimal=round(la, 5), lon_decimal=round(lo, 5), confidence=confidence, note=note))
    return rows, out


def _parse_deg_min(text):
    """'25-31.73N' / '110-59.91W' -> signed decimal degrees."""
    text = str(text).strip()
    hemi = text[-1].upper()
    deg, minutes = text[:-1].split('-')
    value = float(deg) + float(minutes) / 60.0
    return -value if hemi in 'SW' else value


def load_stations(path):
    """Loads ``stations.csv`` (one row per survey) with decimal latitude/longitude.

    Required column: ``survey_id`` (the start of the raw file name, e.g. ``250416WF``).
    Position comes from ``lat_decimal``/``lon_decimal`` if present, otherwise from
    ``lat``/``lon`` written as ``DD-MM.mmN`` / ``DDD-MM.mmW``.  Optional columns
    ``confidence`` and ``note`` are carried through.  Returns an empty frame if the
    file does not exist.
    """
    columns = ['survey_id', 'latitude', 'longitude', 'position_confidence', 'position_note']
    path = pathlib.Path(path)
    if not path.exists():
        return pd.DataFrame(columns=columns)
    df = pd.read_csv(path, dtype=str).fillna('')
    if 'survey_id' not in df:
        raise ValueError(f"{path.name}: missing required column 'survey_id'")
    if {'lat_decimal', 'lon_decimal'} <= set(df.columns):
        lat = pd.to_numeric(df['lat_decimal'], errors='coerce')
        lon = pd.to_numeric(df['lon_decimal'], errors='coerce')
    elif {'lat', 'lon'} <= set(df.columns):
        lat = df['lat'].map(lambda v: _parse_deg_min(v) if v else np.nan)
        lon = df['lon'].map(lambda v: _parse_deg_min(v) if v else np.nan)
    else:
        raise ValueError(f"{path.name}: needs lat_decimal/lon_decimal or lat/lon columns")
    out = pd.DataFrame({
        'survey_id': df['survey_id'].str.strip(),
        'latitude': lat, 'longitude': lon,
        'position_confidence': df.get('confidence', ''),
        'position_note': df.get('note', ''),
    })
    dup = out['survey_id'][out['survey_id'].duplicated()].unique()
    if len(dup):
        raise ValueError(f"{path.name}: duplicate survey_id values: {list(dup)}")
    return out


def survey_id_from_file_key(file_key):
    """'250416WF-D20250416-T191232_1' -> '250416WF' (text before the first '-')."""
    return str(file_key).split('-')[0]


def raw_file_key(path):
    """Filesystem-safe identifier for a raw file, e.g. '250416WF-D20250416-T191232_1'."""
    stem = pathlib.Path(path).stem
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', stem).strip('_')


def ingest_ek80_raw(file_path, nc_dir=None, overwrite=True):
    """Converts one EK80 ``.raw`` file to an Echopype EchoData object.

    The conversion is lossless: the complex samples and all instrument
    settings are repackaged in the SONAR-netCDF4 convention.  When
    ``nc_dir`` is given the converted file is saved there as
    ``<file_key>.nc`` and kept as the archive copy of the raw data.
    """
    file_path = pathlib.Path(file_path)
    echodata = ep.open_raw(file_path, sonar_model='EK80')
    nc_path = None
    if nc_dir is not None:
        nc_dir = pathlib.Path(nc_dir)
        nc_dir.mkdir(parents=True, exist_ok=True)
        nc_path = nc_dir / f"{raw_file_key(file_path)}.nc"
        echodata.to_netcdf(nc_path, overwrite=overwrite)
        logging.info(f"Archived converted data: {nc_path.name}")
    return echodata, nc_path


def sequential_alignment(valid):
    """Pairs up pings when the channels transmit one after another (sequential pinging).

    ``valid`` is a boolean array (channel, ping): True where that channel transmitted.  With
    sequential pinging the file's ping list is the interleaving of every channel's pings, so each
    channel has data only at some of them.  Returns None when every channel pinged every time, or
    (index per channel, common ping times as int64 ns): the k-th ping of each channel is treated as
    simultaneous and given the mean of the channels' times.  The channels' real times differ by
    about one ping interval; that offset is ignored.
    """
    valid = np.asarray(valid, dtype=bool)
    if valid.all():
        return None
    idx = [np.flatnonzero(v) for v in valid]
    n = min(len(i) for i in idx)
    if n == 0:
        raise ValueError("A channel has no pings in this file")
    idx = [i[:n] for i in idx]
    return idx


def _align_dataset(ds, var, idx_and_times=None):
    """Applies :func:`sequential_alignment` to a Dataset with channel and ping_time dimensions."""
    other = [d for d in ds[var].dims if d not in ('channel', 'ping_time')]
    valid = np.isfinite(ds[var]).any(dim=other).transpose('channel', 'ping_time').values
    idx = sequential_alignment(valid)
    if idx is None:
        return ds, None
    t = ds['ping_time'].values.astype('datetime64[ns]').astype('int64')
    common = np.mean([t[i] for i in idx], axis=0).astype('int64').astype('datetime64[ns]')
    parts = [ds.isel(channel=[c], ping_time=i).assign_coords(ping_time=common) for c, i in enumerate(idx)]
    out = xr.concat(parts, dim='channel', data_vars='minimal', coords='minimal', compat='override')
    out.attrs = dict(ds.attrs)
    logging.info(f"Sequential pinging: paired {len(common)} pings per channel (was {ds.sizes['ping_time']} interleaved)")
    return out, idx


def _waveform_and_encode_mode(echodata):
    """Chooses Echopype calibration modes from what was actually transmitted."""
    beam = echodata['Sonar/Beam_group1']
    # With sequential pinging a channel has no value (NaN) at the pings the other channel transmitted.
    types = {str(t) for t in np.unique(beam['transmit_type'].values.astype(object)[pd.notna(beam['transmit_type'].values)])}
    if types <= {'CW'}:
        waveform = 'CW'
    elif types and 'CW' not in types:
        waveform = 'BB'           # LFM / FM broadband pulses (pulse compression)
    else:
        raise ValueError(
            f"File mixes CW and broadband pings ({sorted(types)}); "
            "split the file or process the two waveform modes separately."
        )
    encode = 'complex' if 'backscatter_i' in beam else 'power'
    return waveform, encode


def build_env_params(config):
    """Environment overrides from config.  Empty dict = use the values in the .raw file.

    The EK80 software stores whatever the operator typed in; the defaults
    (10 C, 35 PSU, 100 m) are common when nobody updated them.  Supplying
    measured values changes the sound speed (range) and absorption (Sv).
    """
    env = {}
    if 'TEMPERATURE_C' in config:
        env['temperature'] = float(config['TEMPERATURE_C'])
    if 'SALINITY_PSU' in config:
        env['salinity'] = float(config['SALINITY_PSU'])
    if 'PRESSURE_DBAR' in config:
        env['pressure'] = float(config['PRESSURE_DBAR'])
    if 'SOUND_SPEED_MS' in config:
        env['sound_speed'] = float(config['SOUND_SPEED_MS'])
    elif {'temperature', 'salinity', 'pressure'} <= env.keys():
        # Echopype keeps the file's sound speed unless told otherwise, so
        # derive it from the supplied profile values (Mackenzie, 1981).
        env['sound_speed'] = float(calc_sound_speed(
            temperature=env['temperature'], salinity=env['salinity'],
            pressure=env['pressure'], formula_source='Mackenzie'))
    elif env:
        logging.warning(
            "Partial environment override: set TEMPERATURE_C, SALINITY_PSU and "
            "PRESSURE_DBAR together (or SOUND_SPEED_MS) to change the sound speed; "
            "only absorption will reflect the supplied values.")
    return env


def calibrate_sv(echodata, config):
    """Calibrates complex backscatter to Sv (dB re 1 m-1).

    Applies the sonar equation with the transducer gain, equivalent beam
    angle and sa-correction stored in the file, spreading loss (TVG) and
    frequency-dependent sound absorption.  An optional EK80 calibration
    file (``ECS_FILE`` in config.toml, from a sphere calibration) replaces
    the factory values.
    """
    waveform, encode = _waveform_and_encode_mode(echodata)
    kwargs = dict(waveform_mode=waveform, encode_mode=encode)
    env = build_env_params(config)
    if env:
        kwargs['env_params'] = env
    if config.get('ECS_FILE'):
        kwargs['ecs_file'] = str(config['ECS_FILE'])
    logging.info(f"Calibrating Sv: waveform={waveform}, encode={encode}, "
                 f"env_override={env or 'none (values stored in .raw)'}")
    ds_Sv = ep.calibrate.compute_Sv(echodata, **kwargs)
    ds_Sv = ep.consolidate.add_depth(
        ds_Sv, echodata, depth_offset=float(config.get('TRANSDUCER_DEPTH_M', 0.0)))
    ds_Sv, _ = _align_dataset(ds_Sv, 'Sv')
    ds_Sv['frequency_khz'] = (ds_Sv['frequency_nominal'] / 1000.0).astype(float)
    ds_Sv.attrs['waveform_mode'] = waveform
    ds_Sv.attrs['encode_mode'] = encode
    return ds_Sv


def relative_echo_level(echodata, sound_speed, transducer_depth_m=0.0, block=8):
    """Uncalibrated echo level (dB re arbitrary) for a quick look at the data as recorded.

    Per channel: average the complex samples over the transducer sectors, take power,
    average ``block`` consecutive range samples, and convert to dB.  No gain, spreading
    loss, absorption, noise removal or seabed masking is applied, so values are only
    comparable within one channel.  Returns a Dataset with ``level_db`` (channel, ping_time,
    range_block) and a per-channel depth axis ``depth_m`` (channel, range_block), padded
    with NaN where a channel has fewer samples.
    """
    beam = echodata['Sonar/Beam_group1']
    if 'backscatter_i' not in beam:
        raise ValueError('Uncalibrated echo level needs complex (real + imaginary) samples')
    levels, depths, pinged = [], [], []
    for ch in range(beam.sizes['channel']):
        z = beam['backscatter_r'].isel(channel=ch).values + 1j * beam['backscatter_i'].isel(channel=ch).values
        with np.errstate(all='ignore'):
            power = np.abs(np.nanmean(z, axis=-1)) ** 2          # (ping, range)
        n = power.shape[1] // block * block
        with np.errstate(all='ignore'):
            blocks = np.nanmean(power[:, :n].reshape(power.shape[0], -1, block), axis=2)
        valid = np.isfinite(blocks).any(axis=0)
        last = int(np.where(valid)[0].max()) + 1 if valid.any() else 0
        blocks = blocks[:, :last]
        pinged.append(np.isfinite(blocks).any(axis=1))
        dt = float(np.nanmean(beam['sample_interval'].isel(channel=ch).values))
        range_m = sound_speed * dt * np.arange(n).reshape(-1, block).mean(axis=1)[:last] / 2.0
        levels.append(10 * np.log10(np.maximum(blocks, 1e-30)).astype('float32'))
        depths.append((range_m + transducer_depth_m).astype('float32'))
    times = beam['ping_time'].values
    idx = sequential_alignment(np.array(pinged))
    if idx is not None:                                 # sequential pinging: pair up each channel's pings
        t = times.astype('datetime64[ns]').astype('int64')
        times = np.mean([t[i] for i in idx], axis=0).astype('int64').astype('datetime64[ns]')
        levels = [l[i] for l, i in zip(levels, idx)]
    width = max(l.shape[1] for l in levels)
    pad = lambda a, w: np.pad(a, [(0, 0)] * (a.ndim - 1) + [(0, w - a.shape[-1])], constant_values=np.nan)
    return xr.Dataset(
        {'level_db': (('channel', 'ping_time', 'range_block'), np.stack([pad(l, width) for l in levels])),
         'depth_m': (('channel', 'range_block'), np.stack([pad(d, width) for d in depths]))},
        coords={'channel': beam['channel'].astype(str).values, 'ping_time': times},
        attrs={'description': 'Uncalibrated relative echo level, dB re arbitrary; range blocks of '
                              f'{block} samples'})


def compute_ts(echodata, config):
    """Target strength (dB re 1 m2) per sample, plus split-beam angles where the channel has them.

    Returns (ds_TS, angles) where angles maps channel index -> (along, athwart) arrays of shape
    (ping, range) in degrees.  Single-beam channels (200 kHz here) have no angles.
    """
    waveform, encode = _waveform_and_encode_mode(echodata)
    kwargs = dict(waveform_mode=waveform, encode_mode=encode)
    env = build_env_params(config)
    if env:
        kwargs['env_params'] = env
    if config.get('ECS_FILE'):
        kwargs['ecs_file'] = str(config['ECS_FILE'])
    ds_ts = ep.calibrate.compute_TS(echodata, **kwargs)
    angles = {}
    for i in range(ds_ts.sizes['channel']):
        try:
            a = ep.consolidate.add_splitbeam_angle(ds_ts.isel(channel=[i]).copy(), echodata, waveform_mode=waveform,
                                                   encode_mode=encode, pulse_compression=False, to_disk=False)
            angles[i] = (a['angle_alongship'].isel(channel=0).values, a['angle_athwartship'].isel(channel=0).values)
        except Exception as exc:       # single-beam channels raise "beam_type not recognized"
            logging.debug(f"No split-beam angles for channel {i}: {exc}")
    ds_ts, idx = _align_dataset(ds_ts, 'TS')
    if idx is not None:
        angles = {i: (a[0][idx[i]], a[1][idx[i]]) for i, a in angles.items()}
    return ds_ts, angles


def _per_channel(ds, name):
    """One value per channel for a variable that may or may not vary with ping."""
    if name not in ds:
        return np.full(ds.sizes['channel'], np.nan)
    da = ds[name]
    for dim in [d for d in da.dims if d != 'channel']:
        da = da.mean(dim)
    values = np.asarray(da.values, dtype=float)
    return np.broadcast_to(values, (ds.sizes['channel'],)).copy() if values.ndim == 0 else values


def calibration_table(ds_Sv, echodata, config, cruise_id, file_key):
    """What was used to turn the recording into Sv: water properties, gains and beam settings."""
    beam = echodata['Sonar/Beam_group1']
    env_keys = {'TEMPERATURE_C', 'SALINITY_PSU', 'PRESSURE_DBAR', 'SOUND_SPEED_MS'}
    source = 'config.toml' if env_keys & set(config) else 'values stored in the .raw file'
    dt = np.asarray(beam['sample_interval'].mean('ping_time', skipna=True).values, dtype=float)
    c = _per_channel(ds_Sv, 'sound_speed')
    rows = []
    for i in range(ds_Sv.sizes['channel']):
        rows.append(dict(
            cruise_id=cruise_id, file_key=file_key, channel=str(ds_Sv['channel'].values[i]),
            frequency_khz=float(ds_Sv['frequency_khz'].values[i]), environment_source=source,
            sound_speed_ms=c[i], absorption_db_per_m=_per_channel(ds_Sv, 'sound_absorption')[i],
            temperature_c=_per_channel(ds_Sv, 'temperature')[i], salinity_psu=_per_channel(ds_Sv, 'salinity')[i],
            pressure_dbar=_per_channel(ds_Sv, 'pressure')[i], ph=_per_channel(ds_Sv, 'pH')[i],
            gain_db=_per_channel(ds_Sv, 'gain_correction')[i], sa_correction_db=_per_channel(ds_Sv, 'sa_correction')[i],
            equivalent_beam_angle_db=_per_channel(ds_Sv, 'equivalent_beam_angle')[i],
            beamwidth_along_deg=_per_channel(ds_Sv, 'beamwidth_alongship')[i],
            beamwidth_athwart_deg=_per_channel(ds_Sv, 'beamwidth_athwartship')[i],
            transmit_power_w=float(beam['transmit_power'].isel(channel=i).mean('ping_time', skipna=True).values),
            pulse_duration_s=float(beam['transmit_duration_nominal'].isel(channel=i).mean('ping_time', skipna=True).values),
            tau_effective_s=_per_channel(ds_Sv, 'tau_effective')[i], sample_interval_s=dt[i],
            range_bin_m=c[i] * dt[i] / 2.0,
            impedance_transducer_ohm=_per_channel(ds_Sv, 'impedance_transducer')[i],
            impedance_transceiver_ohm=_per_channel(ds_Sv, 'impedance_transceiver')[i],
            waveform_mode=ds_Sv.attrs.get('waveform_mode', ''), ecs_file=str(config.get('ECS_FILE', '')),
            transducer_depth_m=float(config.get('TRANSDUCER_DEPTH_M', 0.0))))
    return pd.DataFrame(rows)
