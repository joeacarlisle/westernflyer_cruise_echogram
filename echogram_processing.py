"""Echogram processing: noise removal, seabed detection, masking and gridding.

Follows the standard active-acoustics workflow for calibrated Sv data:

    near-field exclusion  ->  background-noise removal (De Robertis & Higginbottom 2007)
    ->  impulse-noise mask (Ryan et al. 2015)  ->  seabed detection
    ->  mask below the seabed  ->  MVBS grid (mean volume backscattering strength)

Every step works on one channel at a time, because the EK80 channels on this
transducer are sampled at different rates (e.g. 3 cm vs 0.6 cm range bins).
Key interface: process_sv() and compute_mvbs().
"""

import logging
import warnings

import echopype as ep
import numpy as np
import pandas as pd
import xarray as xr
from echopype.clean import mask_impulse_noise, remove_background_noise
from echopype.commongrid import compute_MVBS
from echopype.mask import detect_seafloor


def validate_config(config):
    """Fails early on settings that would silently produce wrong results."""
    positive = ['NEAR_FIELD_EXCLUSION_M', 'BN_PING_NUM', 'BN_RANGE_M', 'MVBS_RANGE_BIN_M',
                'SEABED_MIN_RANGE_M']
    for key in positive:
        if key not in config or float(config[key]) <= 0:
            raise ValueError(f"config.toml: '{key}' is required and must be > 0")
    for key in ('MVBS_PING_TIME_BIN',):
        if key not in config:
            raise ValueError(f"config.toml: '{key}' is required (e.g. \"10s\")")
    pd.Timedelta(config['MVBS_PING_TIME_BIN'])   # raises if unparseable
    if 'SEABED_CHANNEL_KHZ' in config and float(config['SEABED_CHANNEL_KHZ']) <= 0:
        raise ValueError("config.toml: 'SEABED_CHANNEL_KHZ' must be > 0 (or leave it out to use each file's lowest frequency)")


def _channel_range_info(ds_ch):
    """Number of valid range samples and the median range-bin size (m) for a 1-channel dataset."""
    er = ds_ch['echo_range'].isel(channel=0, ping_time=0).values
    n_valid = int(np.isfinite(er).sum())
    dr = float(np.nanmedian(np.diff(er[:n_valid])))
    return n_valid, dr


def apply_near_field_exclusion(ds_Sv, config):
    """Masks samples closer to the transducer than NEAR_FIELD_EXCLUSION_M (transmit pulse, near field)."""
    limit = float(config['NEAR_FIELD_EXCLUSION_M'])
    out = ds_Sv.copy()
    out['Sv'] = ds_Sv['Sv'].where(ds_Sv['echo_range'] >= limit)
    return out


def apply_background_noise_removal(ds_Sv, config):
    """Estimates and subtracts background noise per channel; drops samples below the SNR threshold."""
    if not config.get('BACKGROUND_NOISE', True):
        return ds_Sv
    parts = []
    for i in range(ds_Sv.sizes['channel']):
        ds_ch = ds_Sv.isel(channel=[i])
        n_valid, dr = _channel_range_info(ds_ch)
        ds_ch = ds_ch.isel(range_sample=slice(0, n_valid))
        range_sample_num = max(1, int(round(float(config['BN_RANGE_M']) / dr)))
        ping_num = int(config['BN_PING_NUM'])
        logging.debug(f"Background noise ch{i}: ping_num={ping_num}, "
                      f"range_sample_num={range_sample_num} (dr={dr:.4f} m)")
        cleaned = remove_background_noise(
            ds_ch, ping_num=ping_num, range_sample_num=range_sample_num,
            SNR_threshold=f"{float(config.get('SNR_THRESHOLD_DB', 3.0))}dB")
        ds_ch = ds_ch.copy()
        ds_ch['Sv'] = cleaned['Sv_corrected']   # Sv minus noise; NaN where SNR is below threshold
        ds_ch['Sv_noise'] = cleaned['Sv_noise']  # kept for single-target detection and quality control
        parts.append(ds_ch)
    return _concat_channels(parts, ds_Sv)


def apply_impulse_noise_mask(ds_Sv, config):
    """Removes isolated spikes (e.g. interference from other instruments)."""
    if not config.get('IMPULSE_NOISE', True):
        return ds_Sv
    parts = []
    for i in range(ds_Sv.sizes['channel']):
        ds_ch = ds_Sv.isel(channel=[i])
        n_valid, _ = _channel_range_info(ds_ch)
        ds_ch = ds_ch.isel(range_sample=slice(0, n_valid))
        noise = mask_impulse_noise(
            ds_ch, depth_bin=f"{float(config.get('IMPULSE_DEPTH_BIN_M', 5.0))}m",
            num_side_pings=int(config.get('IMPULSE_SIDE_PINGS', 2)),
            impulse_noise_threshold=f"{float(config.get('IMPULSE_THRESHOLD_DB', 10.0))}dB",
            range_var='depth')
        # True in the echopype noise mask means "impulse noise"; apply_mask keeps True.
        keep = ~noise.fillna(0).astype(bool)
        logging.debug(f"Impulse noise ch{i}: flagged {100 * float(1 - keep.mean()):.4f}% of samples")
        keep = keep.transpose(*ds_ch['Sv'].dims)
        if keep.shape != ds_ch['Sv'].shape:
            raise ValueError(f"Impulse-noise mask shape {keep.shape} does not match Sv {ds_ch['Sv'].shape}")
        ds_ch['Sv'] = ds_ch['Sv'].where(keep.values)
        parts.append(ds_ch)
    return _concat_channels(parts, ds_Sv)


def apply_transient_noise_mask(ds_Sv, config):
    """Optional: removes short bursts of noise that are much stronger than the pings around them."""
    if not config.get('TRANSIENT_NOISE', False):
        return ds_Sv
    import echogram_analysis as analysis
    out = ds_Sv.copy()
    sv = out['Sv'].values.copy()                                         # (channel, ping, range)
    for i in range(ds_Sv.sizes['channel']):
        depth = out['depth'].isel(channel=i, ping_time=0).values
        valid = np.isfinite(depth)
        grid = sv[i][:, valid].T                                         # (range, ping)
        noise = analysis.transient_noise_mask(
            grid, depth[valid], depth_bin_m=float(config.get('TRANSIENT_DEPTH_BIN_M', 5.0)),
            side_pings=int(config.get('TRANSIENT_SIDE_PINGS', 25)),
            exclude_above_m=float(config.get('TRANSIENT_EXCLUDE_ABOVE_M', 10.0)),
            threshold_db=float(config.get('TRANSIENT_THRESHOLD_DB', 12.0)))
        logging.debug(f"Transient noise ch{i}: flagged {100 * noise.mean():.4f}% of samples")
        block = sv[i][:, valid]
        block[noise.T] = np.nan
        sv[i][:, valid] = block
    out['Sv'] = (out['Sv'].dims, sv)
    return out


def exclusion_mask(ping_times, config):
    """Boolean per ping: True for pings inside any EXCLUDE_TIME_RANGES_UTC [start, end] pair."""
    ping = pd.to_datetime(np.asarray(ping_times))
    bad = np.zeros(len(ping), bool)
    for pair in config.get('EXCLUDE_TIME_RANGES_UTC', []) or []:
        a, b = pd.Timestamp(pair[0]), pd.Timestamp(pair[1])
        bad |= (ping >= a) & (ping <= b)
    return bad


def apply_exclusions(ds_Sv, config):
    """Blanks pings the user marked as bad (off-transect, interference, ship manoeuvres)."""
    bad = exclusion_mask(ds_Sv['ping_time'].values, config)
    out = ds_Sv.copy()
    if bad.any():
        keep = xr.DataArray(~bad, dims='ping_time', coords={'ping_time': ds_Sv['ping_time']})
        out['Sv'] = ds_Sv['Sv'].where(keep)
        logging.info(f"Excluded {int(bad.sum())} ping(s) listed in EXCLUDE_TIME_RANGES_UTC")
    out.attrs['n_excluded_pings'] = int(bad.sum())
    return out


def _concat_channels(parts, template):
    """Rebuilds the multi-channel dataset (channels may differ in range-sample count)."""
    out = xr.concat(parts, dim='channel', join='outer', data_vars='all')
    out.attrs.update(template.attrs)
    return out


def detect_seabed(ds_Sv, config):
    """Detects the seabed on the configured channel and returns depth (m) per ping.

    Uses Echopype's threshold ('basic') detector: the first sample below
    SEABED_MIN_RANGE_M whose Sv exceeds SEABED_THRESHOLD_DB, minus SEABED_OFFSET_M.
    """
    target = float(config['SEABED_CHANNEL_KHZ'])
    idx = int(np.argmin(np.abs(ds_Sv['frequency_khz'].values - target)))
    ds_ch = ds_Sv.isel(channel=[idx])
    n_valid, dr = _channel_range_info(ds_ch)
    skip_bins = int(np.ceil(float(config['SEABED_MIN_RANGE_M']) / dr))
    channel_name = str(ds_Sv['channel'].values[idx])
    seabed = detect_seafloor(
        ds_ch.isel(range_sample=slice(0, n_valid)), method='basic',
        params=dict(var_name='Sv', channel=channel_name,
                    threshold=float(config.get('SEABED_THRESHOLD_DB', -40.0)),
                    offset_m=float(config.get('SEABED_OFFSET_M', 0.5)),
                    bin_skip_from_surface=skip_bins))
    seabed = seabed.rename('seabed_depth_m')
    seabed.attrs.update(units='m', long_name='Detected seabed depth',
                        detection_frequency_khz=float(ds_Sv['frequency_khz'].values[idx]))
    n_found = int(np.isfinite(seabed).sum())
    logging.info(f"Seabed detected on {seabed.attrs['detection_frequency_khz']:g} kHz: "
                 f"{n_found}/{seabed.size} pings, "
                 f"{float(seabed.min()):.1f}-{float(seabed.max()):.1f} m")
    return seabed


def apply_seabed_mask(ds_Sv, seabed, config):
    """Masks everything at or below the seabed (minus a safety margin) on every channel."""
    margin = float(config.get('SEABED_MASK_MARGIN_M', 1.0))
    depth = ds_Sv['depth']
    keep = depth <= (seabed - margin)
    # Pings with no detection keep all data (comparison with NaN is False, so restore them).
    keep = keep | seabed.isnull()
    out = ds_Sv.copy()
    out['Sv'] = ds_Sv['Sv'].where(keep)
    return out


def process_sv(ds_Sv, config):
    """Runs the full cleaning chain and returns (clean Sv dataset, seabed DataArray)."""
    ds = apply_exclusions(ds_Sv, config)
    ds = apply_near_field_exclusion(ds, config)
    seabed = detect_seabed(ds, config)           # detect on near-field-cleaned, un-noise-subtracted Sv
    bad = exclusion_mask(ds_Sv['ping_time'].values, config)
    if bad.any():                                  # no seabed for pings the user excluded
        attrs = dict(seabed.attrs)
        seabed = seabed.where(xr.DataArray(~bad, dims='ping_time', coords={'ping_time': seabed['ping_time']}))
        seabed.attrs.update(attrs)
    ds = apply_background_noise_removal(ds, config)
    ds = apply_impulse_noise_mask(ds, config)
    ds = apply_transient_noise_mask(ds, config)
    ds = apply_seabed_mask(ds, seabed, config)
    ds['seabed_depth_m'] = seabed
    return ds, seabed


def compute_mvbs(ds_clean, config):
    """Grids cleaned Sv to mean volume backscattering strength (MVBS, dB) cells.

    Averaging is done in the linear domain by Echopype (not on dB values).
    """
    ds_grid = ds_clean[['Sv', 'echo_range', 'depth', 'frequency_nominal']].copy()
    with warnings.catch_warnings():
        # Range padding of the coarser channel is NaN in 'depth'; those cells hold no data, so
        # Echopype's "depth contains NaN" warning is expected here and harmless.
        warnings.filterwarnings('ignore', message='.*depth.*coordinate array contain NaNs.*')
        mvbs = compute_MVBS(
            ds_grid, range_var='depth',
            range_bin=f"{float(config['MVBS_RANGE_BIN_M'])}m",
            ping_time_bin=str(config['MVBS_PING_TIME_BIN']))
    return mvbs


def mvbs_to_frame(mvbs, ds_Sv, cruise_id, file_key):
    """Flattens an MVBS dataset to a long table (one row per non-empty cell) for Parquet/DuckDB."""
    freq = pd.Series(ds_Sv['frequency_khz'].values, index=ds_Sv['channel'].values)
    da = mvbs['Sv'].rename('mvbs_db')
    df = da.to_dataframe().reset_index()
    depth_col = [c for c in df.columns if c in ('depth', 'echo_range')][0]
    df = df.rename(columns={depth_col: 'depth_m'})
    df['channel'] = df['channel'].astype(str)
    df['frequency_khz'] = df['channel'].map(freq).astype(float)
    df = df.dropna(subset=['mvbs_db'])
    df['cruise_id'] = cruise_id
    df['file_key'] = file_key
    return df[['cruise_id', 'file_key', 'channel', 'frequency_khz', 'ping_time', 'depth_m', 'mvbs_db']]


def seabed_to_frame(seabed, cruise_id, file_key):
    """Per-ping seabed depth table."""
    df = seabed.to_dataframe().reset_index()[['ping_time', 'seabed_depth_m']]
    df['frequency_khz'] = float(seabed.attrs.get('detection_frequency_khz', np.nan))
    df['cruise_id'] = cruise_id
    df['file_key'] = file_key
    return df[['cruise_id', 'file_key', 'frequency_khz', 'ping_time', 'seabed_depth_m']]


def qc_summary(ds_cal, ds_clean, seabed, config, cruise_id, file_key):
    """One row per frequency: ping timing, excluded pings, share of data left, noise level."""
    ping = pd.to_datetime(ds_cal['ping_time'].values)
    gaps = np.diff(ping.to_numpy().astype('datetime64[ns]')) / np.timedelta64(1, 's')
    med = float(np.median(gaps)) if len(gaps) else float('nan')
    n_gaps = int((gaps > 1.5 * med).sum()) if len(gaps) else 0
    excluded = exclusion_mask(ds_cal['ping_time'].values, config)
    floor = np.where(np.isfinite(seabed.values), seabed.values, np.inf)[:, None]      # no detection: whole range
    rows = []
    for i in range(ds_clean.sizes['channel']):
        depth = ds_clean['depth'].isel(channel=i).values
        sv = ds_clean['Sv'].isel(channel=i).values
        above = np.isfinite(depth) & (depth <= floor) & (~excluded)[:, None]
        valid_fraction = float(np.isfinite(sv[above]).mean()) if above.any() else float('nan')
        noise = ds_clean['Sv_noise'].isel(channel=i).values if 'Sv_noise' in ds_clean else None
        band = np.isfinite(depth) & (depth >= 40) & (depth <= 60) & (~excluded)[:, None]
        med_noise = float(np.nanmedian(noise[band])) if noise is not None and band.any() else float('nan')
        rows.append(dict(cruise_id=cruise_id, file_key=file_key, frequency_khz=float(ds_clean['frequency_khz'].values[i]),
                         n_pings=int(len(ping)), median_ping_interval_s=med, max_ping_gap_s=float(gaps.max()) if len(gaps) else float('nan'),
                         n_ping_gaps=n_gaps, n_excluded_pings=int(excluded.sum()),
                         valid_fraction_above_seabed=valid_fraction, median_noise_sv_db_40_60m=med_noise))
    return pd.DataFrame(rows)


def find_schools(ds_clean, config, cruise_id, file_key):
    """Detects schools / scattering layers on the block-averaged cleaned Sv of the chosen frequencies."""
    import echogram_analysis as analysis
    wanted = config.get('SCHOOL_FREQUENCIES_KHZ', [38])
    if not config.get('SCHOOL_DETECTION', True) or not wanted:
        return pd.DataFrame()
    frames = []
    for i in range(ds_clean.sizes['channel']):
        f = float(ds_clean['frequency_khz'].values[i])
        if not any(abs(f - float(w)) < 0.5 for w in wanted):
            continue
        depth = ds_clean['depth'].isel(channel=i, ping_time=0).values
        valid = np.isfinite(depth)
        sv = ds_clean['Sv'].isel(channel=i).values[:, valid].T          # (range, ping)
        times, dep, grid = analysis.block_average_sv(
            ds_clean['ping_time'].values, depth[valid], sv,
            float(config.get('SCHOOL_SMOOTH_RANGE_M', 0.5)), int(config.get('SCHOOL_SMOOTH_PINGS', 5)))
        found = analysis.detect_schools(
            times, dep, grid, float(config.get('SCHOOL_THRESHOLD_DB', -70.0)),
            min_thickness_m=float(config.get('SCHOOL_MIN_THICKNESS_M', 2.0)),
            min_length_s=float(config.get('SCHOOL_MIN_LENGTH_S', 5.0)))
        if len(found):
            found.insert(0, 'frequency_khz', f)
            found.insert(0, 'file_key', file_key)
            found.insert(0, 'cruise_id', cruise_id)
            frames.append(found)
    logging.info(f"School/layer detection: {sum(len(x) for x in frames)} detection(s)")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
