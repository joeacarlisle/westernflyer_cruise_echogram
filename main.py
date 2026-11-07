"""Western Flyer EK80 echogram build orchestrator.

Usage:
    python main.py <cruise_id> [--raw-dir PATH] [--db PATH]
                               [--range-bin 1.0] [--ping-bin 10s]
                               [--no-sv-store] [--verbose]

The cruise_id maps to ``cruises/<cruise_id>/`` and is the authoritative
source of truth for the cruise.  Expected layout:

    cruises/<cruise_id>/
        config.toml     processing parameters (environment, noise, seabed, MVBS grid)
        stations.csv    optional: survey_id, lat, lon (the .raw files carry no GPS here)
        nav.csv         optional: time_utc, lat, lon ship log; gives each ping a position, speed and distance
        raw/            raw Simrad EK80 .raw files - one per recording

Processing pipeline (per .raw file):
    ingest_ek80_raw         .raw -> EchoData (lossless SONAR-netCDF4 archive .nc only with --write-netcdf)
    calibrate_sv            complex backscatter -> Sv (dB re 1 m-1), add depth
    process_sv              excluded pings, near-field exclusion, seabed detection, background-noise
                            removal, impulse-noise mask, optional transient-noise mask, mask below seabed
    compute_mvbs            grid cleaned Sv to MVBS cells (range bin x ping-time bin)
    compute_ts + detect_single_targets   target strength, split-beam angles, single-target detection
    find_schools            school / scattering-layer detection on the block-averaged cleaned Sv
    calibration_table, qc_summary        what was used for calibration; ping timing and data quality
    interpolate_track       optional ship track from cruises/<cruise_id>/nav.csv

Outputs (all under processed/):
    processed/<cruise_id>/netcdf/<file_key>.nc     converted raw data (archive; only with --write-netcdf)
    processed/<cruise_id>/sv/<file_key>.zarr       cleaned full-resolution Sv (dashboard)
    processed/<cruise_id>/raw_level/<file_key>.zarr uncalibrated echo level, both frequencies (dashboard)
    processed/<cruise_id>/mvbs.parquet, seabed.parquet, files.parquet   (only with --write-parquet)
    processed/wf_echogram.duckdb
        echogram_mvbs     - one row per (file, channel, ping-time bin, depth bin)
        echogram_seabed   - one row per ping
        echogram_files    - one row per raw file (with survey position from stations.csv)
        echogram_stations - one row per survey in stations.csv, flagged if a raw file matched
        echogram_calibration - water properties, gains and beam settings used, per file and frequency
        echogram_qc       - ping timing, excluded pings, data left after cleaning, noise level
        echogram_targets  - detected single targets (TS, depth, beam position)
        echogram_schools  - detected schools / scattering layers
        echogram_track    - per-ping position, distance and speed (only when nav.csv exists)
        echogram_config   - the config.toml used for each build (for the Processing tab)
        build_metadata    - one row per build, with SHA-256 of config.toml and raw files

Per-cruise records are replaced (not appended) inside a single
transaction, so a failed insert can't leave the database empty.
"""

import argparse
import datetime as dt
import hashlib
import json
import logging
import pathlib
import subprocess
import sys

import duckdb
import echopype
import numpy as np
import pandas as pd
import xarray as xr

import echogram_processing
import ek80_ingestion
import echogram_analysis
import echogram_targets
from ek80_ingestion import (calibrate_sv, calibration_table, compute_ts, ingest_ek80_raw, load_stations,
                            raw_file_key, relative_echo_level, survey_id_from_file_key)

BASE_DIR = pathlib.Path(__file__).resolve().parent
DEFAULT_DB_PATH = BASE_DIR / "processed" / "wf_echogram.duckdb"
LOG_DIR = BASE_DIR / "logs"


def setup_logging(cruise_id, verbose=False):
    """Sets up logging with a per-build, timestamped log file."""
    LOG_DIR.mkdir(exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = LOG_DIR / f"wf_echogram_build_{cruise_id}_{timestamp}.log"
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path, mode='w'),
            logging.StreamHandler(sys.stdout),
        ],
    )
    # Third-party libraries are chatty at DEBUG; keep their output at WARNING.
    for noisy in ('numcodecs', 'flox', 'matplotlib', 'fsspec', 'asyncio'):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    return log_path


def file_sha256(path):
    """Returns the hex SHA-256 of a file (used for provenance)."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def git_commit():
    """Returns the current git commit hash, or 'unknown' if not a repo."""
    try:
        out = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=BASE_DIR, stderr=subprocess.DEVNULL)
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 'unknown'


def write_sv_store(ds_clean, zarr_path):
    """Writes the cleaned full-resolution Sv as a chunked Zarr store for lazy dashboard reads.

    Only what the dashboard needs is kept: Sv (float32, NaN = removed), a per-channel
    display depth axis, and the seabed line.
    """
    ds_clean = ds_clean.assign_coords(channel=ds_clean['channel'].astype(str))
    display_depth = ds_clean['depth'].isel(ping_time=0).drop_vars('ping_time', errors='ignore')
    store = xr.Dataset(
        {
            'Sv': ds_clean['Sv'].astype('float32'),
            'display_depth_m': display_depth.astype('float32'),
            'seabed_depth_m': ds_clean['seabed_depth_m'].astype('float32'),
            'frequency_khz': ds_clean['frequency_khz'].astype('float32'),
        },
        attrs={'description': 'Cleaned Sv (dB re 1 m-1); NaN = removed or below seabed',
               'echopype_version': echopype.__version__},
    )
    store['Sv'] = store['Sv'].chunk({'channel': 1, 'ping_time': 128, 'range_sample': 4096})
    store.to_zarr(zarr_path, mode='w', consolidated=False)


# Column types of the tables, so that an empty result (for example no targets detected) still creates a
# table with the right types and later builds with data do not clash with it.
TABLE_TYPES = {
    'echogram_calibration': {'cruise_id': 'str', 'file_key': 'str', 'channel': 'str', 'environment_source': 'str',
                             'waveform_mode': 'str', 'ecs_file': 'str'},
    'echogram_qc': {'cruise_id': 'str', 'file_key': 'str', 'n_pings': 'int64', 'n_ping_gaps': 'int64',
                    'n_excluded_pings': 'int64'},
    'echogram_targets': {'cruise_id': 'str', 'file_key': 'str', 'ping_index': 'int64', 'ping_time': 'datetime64[ns]',
                         'compensated': 'bool'},
    'echogram_schools': {'cruise_id': 'str', 'file_key': 'str', 'school_id': 'int64', 'start_time': 'datetime64[ns]',
                         'end_time': 'datetime64[ns]', 'n_cells': 'int64'},
    'echogram_track': {'cruise_id': 'str', 'file_key': 'str', 'ping_time': 'datetime64[ns]'},
}
TABLE_COLUMNS = {
    'echogram_targets': ['cruise_id', 'file_key'] + echogram_targets.TARGET_COLUMNS,
    'echogram_schools': ['cruise_id', 'file_key', 'frequency_khz', 'school_id', 'start_time', 'end_time', 'duration_s',
                         'depth_top_m', 'depth_bottom_m', 'thickness_m', 'mean_depth_m', 'mean_sv_db', 'max_sv_db',
                         'n_cells'],
    'echogram_track': ['cruise_id', 'file_key', 'ping_time', 'lat', 'lon', 'distance_m', 'speed_kn'],
}


def typed_table(name, df):
    """Gives a result frame fixed column names and types (floats unless listed in TABLE_TYPES)."""
    cols = TABLE_COLUMNS.get(name)
    if cols is not None:
        df = df.reindex(columns=cols)
    types = TABLE_TYPES.get(name, {})
    for col in df.columns:
        want = types.get(col, 'float64')
        if want == 'str':
            df[col] = df[col].astype('object').where(df[col].notna(), '').astype(str)
        elif want == 'bool':
            df[col] = df[col].fillna(False).astype(bool)
        elif want == 'int64':
            df[col] = df[col].fillna(0).astype('int64')
        elif want.startswith('datetime'):
            df[col] = pd.to_datetime(df[col])
        else:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('float64')
    return df.reset_index(drop=True)


def detect_targets_for_file(echodata, ds_Sv, ds_clean, seabed, config, cruise_id, file_key):
    """Computes TS (and split-beam angles) and runs the single-target detector on every channel."""
    if not config.get('TARGET_DETECTION', True):
        return pd.DataFrame()
    ds_ts, angles = compute_ts(echodata, config)
    beam = echodata['Sonar/Beam_group1']
    margin = float(config.get('TARGET_SEABED_MARGIN_M', 2.0))
    max_depth = (seabed.values - margin) if seabed.notnull().all() else None
    bad_pings = echogram_processing.exclusion_mask(ds_Sv['ping_time'].values, config)
    frames = []
    for i in range(ds_Sv.sizes['channel']):
        depth = ds_Sv['depth'].isel(channel=i, ping_time=0).values
        dt_s = float(np.nanmean(beam['sample_interval'].isel(channel=i).values))   # NaN at pings the other channel transmitted
        tau = float(np.nanmean(ds_Sv['tau_effective'].isel(channel=i).values)) if 'tau_effective' in ds_Sv else \
            float(np.nanmean(beam['transmit_duration_nominal'].isel(channel=i).values))
        along, athwart = angles.get(i, (None, None))
        bw_a = float(np.nanmean(ds_Sv['beamwidth_alongship'].isel(channel=i).values)) if 'beamwidth_alongship' in ds_Sv else None
        bw_b = float(np.nanmean(ds_Sv['beamwidth_athwartship'].isel(channel=i).values)) if 'beamwidth_athwartship' in ds_Sv else None
        noise = ds_clean['Sv_noise'].isel(channel=i).values
        found = echogram_targets.detect_single_targets(
            ds_ts['TS'].isel(channel=i).values, ds_Sv['Sv'].isel(channel=i).values, noise, depth,
            ds_Sv['ping_time'].values, tau / dt_s, float(ds_Sv['frequency_khz'].values[i]),
            along=along, athwart=athwart, bw_along=bw_a, bw_athwart=bw_b,
            ts_min_db=float(config.get('TARGET_TS_MIN_DB', -60.0)), snr_min_db=float(config.get('TARGET_SNR_MIN_DB', 10.0)),
            min_length=float(config.get('TARGET_MIN_LENGTH', 0.8)), max_length=float(config.get('TARGET_MAX_LENGTH', 1.8)),
            max_beam_comp_db=float(config.get('TARGET_MAX_BEAM_COMP_DB', 6.0)),
            min_depth_m=float(config.get('TARGET_MIN_DEPTH_M', 5.0)), max_depth_m=max_depth,
            background_window_m=float(config.get('TARGET_BACKGROUND_WINDOW_M', 5.0)),
            min_peak_over_background_db=float(config.get('TARGET_PEAK_OVER_BACKGROUND_DB', 15.0)))
        if bad_pings.any() and len(found):
            found = found[~bad_pings[found['ping_index'].to_numpy()]]
        frames.append(found)
    out = pd.concat([f for f in frames if len(f)], ignore_index=True) if any(len(f) for f in frames) else pd.DataFrame()
    logging.info(f"Single targets detected: {len(out)}")
    if len(out):
        out.insert(0, 'file_key', file_key)
        out.insert(0, 'cruise_id', cruise_id)
    return out


def process_one_file(file_path, config, cruise_id, cruise_out, nav=None, write_sv=True, write_nc=False):
    """Runs the full per-file pipeline.  Returns (tables, files_row); tables maps table name -> DataFrame."""
    file_key = raw_file_key(file_path)
    echodata, nc_path = ingest_ek80_raw(file_path, nc_dir=cruise_out / 'netcdf' if write_nc else None)
    ds_Sv = calibrate_sv(echodata, config)
    config = ek80_ingestion.resolve_channel_settings(config, ds_Sv)     # channel choices left out of the config
    ds_clean, seabed = echogram_processing.process_sv(ds_Sv, config)

    mvbs = echogram_processing.compute_mvbs(ds_clean, config)
    tables = {
        'echogram_mvbs': echogram_processing.mvbs_to_frame(mvbs, ds_Sv, cruise_id, file_key),
        'echogram_seabed': echogram_processing.seabed_to_frame(seabed, cruise_id, file_key),
        'echogram_calibration': calibration_table(ds_Sv, echodata, config, cruise_id, file_key),
        'echogram_qc': echogram_processing.qc_summary(ds_Sv, ds_clean, seabed, config, cruise_id, file_key),
        'echogram_schools': echogram_processing.find_schools(ds_clean, config, cruise_id, file_key),
        'echogram_targets': detect_targets_for_file(echodata, ds_Sv, ds_clean, seabed, config, cruise_id, file_key),
    }
    if nav is not None and len(nav):
        track = echogram_analysis.interpolate_track(nav, ds_Sv['ping_time'].values)
        if track['lat'].notna().any():
            track.insert(0, 'file_key', file_key)
            track.insert(0, 'cruise_id', cruise_id)
            tables['echogram_track'] = track
            logging.info(f"Track: {int(track['lat'].notna().sum())}/{len(track)} pings positioned from nav.csv")
        else:
            logging.warning(f"nav.csv does not cover {file_path.name}; no track for this file")
    tables = {k: (typed_table(k, v) if k in TABLE_TYPES else v) for k, v in tables.items()}

    zarr_path = None
    if write_sv:
        zarr_path = cruise_out / 'sv' / f"{file_key}.zarr"
        zarr_path.parent.mkdir(parents=True, exist_ok=True)
        write_sv_store(ds_clean, zarr_path)
        logging.info(f"Wrote cleaned Sv: {zarr_path.name}")

    raw_level_path = None
    if write_sv:
        try:
            level = relative_echo_level(
                echodata, sound_speed=float(np.nanmean(ds_Sv['sound_speed'].values)),
                transducer_depth_m=float(config.get('TRANSDUCER_DEPTH_M', 0.0)))
            raw_level_path = cruise_out / 'raw_level' / f"{file_key}.zarr"
            raw_level_path.parent.mkdir(parents=True, exist_ok=True)
            level.to_zarr(raw_level_path, mode='w', consolidated=False)
            logging.info(f"Wrote uncalibrated echo level: {raw_level_path.name}")
        except Exception as exc:     # the preview is optional; never fail the build for it
            logging.warning(f"Skipped uncalibrated echo-level preview for {file_path.name}: {exc}")
            raw_level_path = None

    ping_times = pd.to_datetime(ds_Sv['ping_time'].values)
    files_row = {
        'cruise_id': cruise_id,
        'file_key': file_key,
        'filename': file_path.name,
        'raw_sha256': file_sha256(file_path),
        'raw_bytes': file_path.stat().st_size,
        'first_ping': ping_times.min(),
        'last_ping': ping_times.max(),
        'n_pings': int(len(ping_times)),
        'frequencies_khz': ','.join(f"{f:g}" for f in ds_Sv['frequency_khz'].values),
        'waveform_mode': ds_Sv.attrs.get('waveform_mode', ''),
        'sound_speed_ms': float(np.nanmean(ds_Sv['sound_speed'].values)),
        'seabed_min_m': float(np.nanmin(seabed.values)) if seabed.notnull().any() else np.nan,
        'seabed_max_m': float(np.nanmax(seabed.values)) if seabed.notnull().any() else np.nan,
        'netcdf_path': str(nc_path.relative_to(BASE_DIR)) if nc_path else '',
        'sv_zarr_path': str(zarr_path.relative_to(BASE_DIR)) if zarr_path else '',
        'raw_level_zarr_path': str(raw_level_path.relative_to(BASE_DIR)) if raw_level_path else '',
    }
    return tables, files_row


def _type_position_columns(df):
    """Position columns get fixed types even when no position is known.

    With no stations.csv every value is missing, and an untyped all-missing column becomes DuckDB's
    NULL type, which later cannot be filled with real positions.
    """
    for col in ('latitude', 'longitude'):
        df[col] = pd.to_numeric(df[col], errors='coerce').astype('float64')
    for col in ('position_confidence', 'position_note'):
        df[col] = df[col].astype('object').where(df[col].notna(), '').astype(str)
    return df


def attach_positions(files_df, stations_df, cruise_id):
    """Adds survey_id, latitude, longitude and position metadata to the file table.

    Returns (files_df, stations_out) where stations_out lists every survey in
    stations.csv with a ``processed`` flag.  Files without a matching survey get
    NaN positions and a warning.
    """
    files_df = files_df.copy()
    files_df['survey_id'] = files_df['file_key'].map(survey_id_from_file_key)
    pos_cols = ['survey_id', 'latitude', 'longitude', 'position_confidence', 'position_note']
    files_df = files_df.merge(stations_df[pos_cols], on='survey_id', how='left')
    files_df = _type_position_columns(files_df)
    for row in files_df[files_df['latitude'].isna()].itertuples():
        logging.warning(f"No position in stations.csv for survey '{row.survey_id}' ({row.filename})")
    stations_out = _type_position_columns(stations_df[pos_cols].copy())
    stations_out['processed'] = stations_out['survey_id'].isin(files_df['survey_id'])
    stations_out.insert(0, 'cruise_id', cruise_id)
    return files_df, stations_out


def _ensure_table(con, table, view):
    """Creates ``table`` from ``view``'s schema; if it exists with different columns, rebuilds it.

    A schema change (new column) would otherwise make the INSERT fail.  Rebuilding drops the
    old rows of *other* cruises in that table, so re-run those builds afterwards.
    """
    exists = con.execute(
        "SELECT count(*) FROM information_schema.tables WHERE table_name = ?", [table]).fetchone()[0]
    if exists:
        old = [r[0] for r in con.execute(f"DESCRIBE {table}").fetchall()]
        new = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {view}").fetchall()]
        if old == new:
            # Same columns: repair any column an earlier build left typed NULL (it was all-missing then).
            old_t = dict(con.execute(f"SELECT column_name, data_type FROM information_schema.columns "
                                     f"WHERE table_name = '{table}'").fetchall())
            for row in con.execute(f"DESCRIBE SELECT * FROM {view}").fetchall():
                col, typ = row[0], row[1]
                if str(old_t.get(col)).strip('"') == 'NULL' and typ.strip('"') != 'NULL':
                    logging.warning(f"Column {table}.{col} was typed NULL by an earlier build; changing it to {typ}.")
                    con.execute(f'ALTER TABLE {table} ALTER COLUMN "{col}" SET DATA TYPE {typ}')
        else:
            logging.warning(f"Table {table} has a different schema; recreating it. "
                            "Rebuild any other cruises to repopulate it.")
            con.execute(f"DROP TABLE {table}")
    con.execute(f"CREATE TABLE IF NOT EXISTS {table} AS SELECT * FROM {view} WHERE 1=0")


def commit_to_duckdb(tables, cruise_id, db_path, config_path, config, build_started):
    """Replaces the cruise's rows in every table atomically and records build provenance.

    ``tables`` maps table name -> DataFrame (all with a cruise_id column).
    """
    db_path.parent.mkdir(exist_ok=True)
    config_hash = file_sha256(config_path)
    files_df = tables['echogram_files']
    manifest = hashlib.sha256(
        ''.join(sorted(files_df['raw_sha256'])).encode()).hexdigest()
    commit = git_commit()
    config_json = json.dumps({k: v for k, v in config.items()}, default=str, sort_keys=True)
    config_df = pd.DataFrame([{'cruise_id': cruise_id, 'config_json': config_json,
                               'config_toml': config_path.read_text(encoding='utf-8')}])
    tables = {**tables, 'echogram_config': config_df}

    with duckdb.connect(str(db_path)) as con:
        for name, df in tables.items():
            con.register(f"{name}_view", df)
            _ensure_table(con, name, f"{name}_view")
        con.execute("""
            CREATE TABLE IF NOT EXISTS build_metadata (
                cruise_id          VARCHAR,
                build_started_ts   TIMESTAMP,
                build_finished_ts  TIMESTAMP,
                config_sha256      VARCHAR,
                raw_manifest_sha256 VARCHAR,
                echopype_version   VARCHAR,
                git_commit         VARCHAR,
                files_processed    BIGINT,
                mvbs_rows_written  BIGINT
            )
        """)
        con.execute("BEGIN TRANSACTION")
        try:
            for name in tables:
                con.execute(f"DELETE FROM {name} WHERE cruise_id = ?", [cruise_id])
                con.execute(f"INSERT INTO {name} SELECT * FROM {name}_view")
            if 'echogram_track' not in tables and con.execute(
                    "SELECT count(*) FROM information_schema.tables WHERE table_name = 'echogram_track'").fetchone()[0]:
                con.execute("DELETE FROM echogram_track WHERE cruise_id = ?", [cruise_id])   # nav.csv was removed
            con.execute("DELETE FROM build_metadata WHERE cruise_id = ?", [cruise_id])
            con.execute(
                "INSERT INTO build_metadata VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [cruise_id, build_started, dt.datetime.now(), config_hash, manifest,
                 echopype.__version__, commit, len(files_df), len(tables['echogram_mvbs'])],
            )
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build the Western Flyer EK80 echogram DuckDB for one cruise.")
    parser.add_argument("cruise_id", help="Cruise directory name under cruises/, e.g. 20250416_Baja")
    parser.add_argument("--range-bin", type=float, default=None,
                        help="Override MVBS_RANGE_BIN_M from config.toml (m)")
    parser.add_argument("--ping-bin", type=str, default=None,
                        help="Override MVBS_PING_TIME_BIN from config.toml, e.g. 10s")
    parser.add_argument("--db", type=pathlib.Path, default=DEFAULT_DB_PATH,
                        help=f"DuckDB output path (default: {DEFAULT_DB_PATH})")
    parser.add_argument("--raw-dir", type=pathlib.Path, default=None,
                        help="Directory containing raw .raw files "
                             "(default: cruises/<cruise_id>/raw/)")
    parser.add_argument("--no-sv-store", action="store_true",
                        help="Skip writing the full-resolution Sv Zarr stores (MVBS only)")
    parser.add_argument("--write-netcdf", action="store_true",
                        help="Also save each .raw file as a lossless SONAR-netCDF4 archive (.nc); off by default")
    parser.add_argument("--write-parquet", action="store_true",
                        help="Also write the tables as Parquet files beside the database; off by default")
    parser.add_argument("--init-config", action="store_true",
                        help="Write a starter cruises/<cruise_id>/config.toml from the first .raw file "
                             "(frequencies found, water values stored in the file) and exit; "
                             "never overwrites an existing file")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable DEBUG-level logging")
    return parser.parse_args(argv)


DEFAULTS_CONFIG_PATH = BASE_DIR / 'config.defaults.toml'


def write_stations_from_gps(cruise_dir, files):
    """No stations.csv: build one from the GPS fixes recorded inside the raw files (never overwrites).

    A survey whose ship moved also gets a time-based ``nav.csv`` (unless one exists) so the dashboard
    can draw the track.  Nothing is written when the files hold no usable GPS.
    """
    try:
        rows, fixes = ek80_ingestion.stations_from_gps(files)
    except Exception as exc:
        logging.warning(f"Could not read GPS from the raw files: {exc}")
        return None
    if not rows:
        logging.warning("No stations.csv, and the raw files hold no usable GPS: positions will be empty.")
        return None
    target = cruise_dir / 'stations.csv'
    cols = ['survey_id', 'lat', 'lon', 'lat_decimal', 'lon_decimal', 'confidence', 'note']
    pd.DataFrame(rows)[cols].to_csv(target, index=False)
    logging.warning(f"No stations.csv: wrote {target} from GPS fixes inside the raw files "
                    f"({', '.join(r['survey_id'] for r in rows)}). Check the positions.")
    moved = [r['survey_id'] for r in rows if 'MOVED' in r['note']]
    nav_path = cruise_dir / 'nav.csv'
    if moved and not nav_path.exists():
        track = pd.concat([f for f in fixes.values()], ignore_index=True).sort_values('time')
        track['time_utc'] = track['time'].dt.floor('s')
        track = track.drop_duplicates('time_utc')[['time_utc', 'lat', 'lon']].round({'lat': 6, 'lon': 6})
        track.to_csv(nav_path, index=False)
        logging.warning(f"The ship moved during {', '.join(moved)}: wrote a ship track to {nav_path}.")
    return target


def _toml_value(v):
    """A Python value as TOML text."""
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (list, tuple)):
        return '[' + ', '.join(_toml_value(x) for x in v) + ']'
    return repr(v)


def write_default_config(cruise_dir, raw_dir, auto=False):
    """Writes ``cruises/<id>/config.toml`` with every setting at its default, from the first .raw file.

    The file lists the channels found, the water values stored in the file and the seabed channel chosen
    (lowest frequency), and carries ``CONFIG_CHECKED = false``: the settings are untested defaults until
    you have reviewed them on this cruise's data and changed that to true.  Never overwrites.
    Returns (path or None, notes list).
    """
    target = cruise_dir / 'config.toml'
    if target.exists():
        return None, [f"{target} already exists; not overwriting it."]
    files = sorted(p for p in raw_dir.iterdir() if p.suffix.lower() == '.raw') if raw_dir.exists() else []
    if not files:
        return None, [f"No .raw files in {raw_dir}; add them first."]
    ed = echopype.open_raw(str(files[0]), sonar_model='EK80')
    beam, env = ed['Sonar/Beam_group1'], ed['Environment']
    freqs = [float(f) / 1000.0 for f in beam['frequency_nominal'].values]
    names = [str(c) for c in beam['channel'].values]
    types = sorted({str(t) for t in np.unique(beam['transmit_type'].values.astype(object)[pd.notna(beam['transmit_type'].values)])})
    sequential = bool(np.isnan(beam['sample_interval'].values.astype(float)).any())

    def first(name):
        try:
            return float(np.ravel(env[name].values)[0])
        except Exception:
            return None
    t, sal, p, c = first('temperature'), first('salinity'), first('depth'), first('sound_speed_indicative')
    lowest = min(freqs)
    notes = []
    if (t, sal, p) == (10.0, 35.0, 100.0):
        notes.append("The water values stored in the file are the EK80 factory defaults (10 C, 35 PSU, 100 m), so they were "
                     "probably never set: measure or look them up and set them in [environment].")
    if len(types) > 1:
        notes.append(f"WARNING: the file mixes pulse types {types}; mixed files are rejected - split them first.")
    if sequential:
        notes.append("Sequential pinging detected: each channel's pings are paired in order and treated as simultaneous.")

    d = ek80_ingestion.load_config(DEFAULTS_CONFIG_PATH) if DEFAULTS_CONFIG_PATH.exists() else {}
    v = lambda key, fallback: _toml_value(d.get(key, fallback))
    out = [f"""# Western Flyer EK80 echogram - cruise configuration, written automatically with DEFAULT SETTINGS
# from {files[0].name}
#
# *** CHECK THESE SETTINGS BEFORE YOU TRUST THE RESULTS. ***
# Every value below is a default, not a value tuned for this cruise's data.  Check, on the echogram
# and in the dashboard's Calibration & QC and Processing tabs:
#   - the red seabed line follows the real seabed (SEABED_CHANNEL_KHZ, SEABED_THRESHOLD_DB)
#   - noisy channels/depths are removed (SNR_THRESHOLD_DB, BN_*), and real echoes are not
#   - the water values (below) are measured ones, not the EK80 software defaults
# When you have checked them, change CONFIG_CHECKED to true; the warnings then stop.
# Keys left out of this file fall back to config.defaults.toml in the project folder.
# Section names are for readability only; every key must be unique across sections.

[metadata]
SONAR = {json.dumps('EK80, channels: ' + '; '.join(f'{n} ({f:g} kHz)' for n, f in zip(names, freqs)))}
CONFIG_CHECKED = false
# Pulse type(s) in the file: {', '.join(types) or 'unknown'}"""]
    out += ['# NOTE: ' + n for n in notes]
    out.append(f"""
[environment]
# Water properties used for sound speed (range) and sound absorption (Sv).  The values stored in the
# file (what was typed into the EK80 software) are used unless you set all four below.
# Stored in the file: temperature {t} C, salinity {sal} PSU, depth {p} m, sound speed {c}.
# TEMPERATURE_C  = {t if t is not None else 10.0}
# SALINITY_PSU   = {sal if sal is not None else 35.0}
# PRESSURE_DBAR  = {p if p is not None else 10.0}
# SOUND_SPEED_MS = {round(c, 1) if c is not None else 1500.0}

[calibration]
# Optional EK80 calibration (.ecs) file from a calibration-sphere run (otherwise the file's factory values are used).
# ECS_FILE = "cruises/{cruise_dir.name}/calibration.ecs"

[platform]
TRANSDUCER_DEPTH_M = {v('TRANSDUCER_DEPTH_M', 0.0)}

[processing]
NEAR_FIELD_EXCLUSION_M = {v('NEAR_FIELD_EXCLUSION_M', 3.0)}
BACKGROUND_NOISE = {v('BACKGROUND_NOISE', True)}
BN_PING_NUM      = {v('BN_PING_NUM', 20)}
BN_RANGE_M       = {v('BN_RANGE_M', 5.0)}
SNR_THRESHOLD_DB = {v('SNR_THRESHOLD_DB', 3.0)}
IMPULSE_NOISE        = {v('IMPULSE_NOISE', True)}
IMPULSE_DEPTH_BIN_M  = {v('IMPULSE_DEPTH_BIN_M', 5.0)}
IMPULSE_SIDE_PINGS   = {v('IMPULSE_SIDE_PINGS', 2)}
IMPULSE_THRESHOLD_DB = {v('IMPULSE_THRESHOLD_DB', 10.0)}
# Seabed: the lowest-frequency channel of the first file ({lowest:g} kHz).
SEABED_CHANNEL_KHZ   = {lowest:g}
SEABED_MIN_RANGE_M   = {v('SEABED_MIN_RANGE_M', 5.0)}
SEABED_THRESHOLD_DB  = {v('SEABED_THRESHOLD_DB', -40.0)}
SEABED_OFFSET_M      = {v('SEABED_OFFSET_M', 0.5)}
SEABED_MASK_MARGIN_M = {v('SEABED_MASK_MARGIN_M', 1.0)}
TRANSIENT_NOISE           = {v('TRANSIENT_NOISE', False)}
TRANSIENT_DEPTH_BIN_M     = {v('TRANSIENT_DEPTH_BIN_M', 5.0)}
TRANSIENT_SIDE_PINGS      = {v('TRANSIENT_SIDE_PINGS', 25)}
TRANSIENT_EXCLUDE_ABOVE_M = {v('TRANSIENT_EXCLUDE_ABOVE_M', 10.0)}
TRANSIENT_THRESHOLD_DB    = {v('TRANSIENT_THRESHOLD_DB', 12.0)}
# Pings inside any [start, end] pair (UTC) are blanked on every channel, e.g. [["2025-04-16T19:13:00", "2025-04-16T19:13:05"]]
EXCLUDE_TIME_RANGES_UTC = []

[targets]
TARGET_DETECTION         = {v('TARGET_DETECTION', True)}
TARGET_TS_MIN_DB         = {v('TARGET_TS_MIN_DB', -60.0)}
TARGET_SNR_MIN_DB        = {v('TARGET_SNR_MIN_DB', 10.0)}
TARGET_MIN_LENGTH        = {v('TARGET_MIN_LENGTH', 0.8)}
TARGET_MAX_LENGTH        = {v('TARGET_MAX_LENGTH', 1.8)}
TARGET_MAX_BEAM_COMP_DB  = {v('TARGET_MAX_BEAM_COMP_DB', 6.0)}
TARGET_BACKGROUND_WINDOW_M      = {v('TARGET_BACKGROUND_WINDOW_M', 5.0)}
TARGET_PEAK_OVER_BACKGROUND_DB  = {v('TARGET_PEAK_OVER_BACKGROUND_DB', 15.0)}
TARGET_MIN_DEPTH_M       = {v('TARGET_MIN_DEPTH_M', 5.0)}
TARGET_SEABED_MARGIN_M   = {v('TARGET_SEABED_MARGIN_M', 2.0)}

[schools]
SCHOOL_DETECTION      = {v('SCHOOL_DETECTION', True)}
SCHOOL_FREQUENCIES_KHZ = [{lowest:g}]
SCHOOL_THRESHOLD_DB   = {v('SCHOOL_THRESHOLD_DB', -70.0)}
SCHOOL_MIN_THICKNESS_M = {v('SCHOOL_MIN_THICKNESS_M', 2.0)}
SCHOOL_MIN_LENGTH_S   = {v('SCHOOL_MIN_LENGTH_S', 5.0)}
SCHOOL_SMOOTH_RANGE_M = {v('SCHOOL_SMOOTH_RANGE_M', 0.5)}
SCHOOL_SMOOTH_PINGS   = {v('SCHOOL_SMOOTH_PINGS', 5)}

[output]
MVBS_RANGE_BIN_M   = {v('MVBS_RANGE_BIN_M', 1.0)}
MVBS_PING_TIME_BIN = {v('MVBS_PING_TIME_BIN', '10s')}
""")
    target.write_text('\n'.join(out), encoding='utf-8')
    return target, notes


def init_config(cruise_dir, raw_dir):
    """--init-config: writes the default-settings config.toml and prints what to check."""
    target, notes = write_default_config(cruise_dir, raw_dir)
    if target is None:
        print(notes[0])
        return 1
    print(f"Wrote {target} (all settings at their defaults, CONFIG_CHECKED = false)")
    for n in notes:
        print("  " + n)
    print("  CHECK the settings on the echogram and in the dashboard's Calibration & QC tab, then set CONFIG_CHECKED = true.")
    return 0


def main(argv=None):
    args = parse_args(argv)
    cruise_id = args.cruise_id
    cruise_dir = BASE_DIR / "cruises" / cruise_id
    if args.init_config:
        sys.exit(init_config(cruise_dir, args.raw_dir if args.raw_dir else cruise_dir / 'raw'))

    (BASE_DIR / "processed").mkdir(exist_ok=True)
    log_path = setup_logging(cruise_id, verbose=args.verbose)

    if not cruise_dir.exists():
        logging.error(f"Cruise directory '{cruise_dir}' not found.")
        sys.exit(1)

    build_started = dt.datetime.now()
    logging.info(f"--- STARTING ECHOGRAM PIPELINE FOR {cruise_id.upper()} ---")
    logging.info(f"Build log: {log_path}")
    logging.info(f"Echopype {echopype.__version__}")

    cruise_config_path = cruise_dir / 'config.toml'
    raw_dir = args.raw_dir if args.raw_dir else cruise_dir / 'raw'
    if not cruise_config_path.exists():
        try:
            written, _ = write_default_config(cruise_dir, raw_dir, auto=True)
        except Exception as exc:
            written = None
            logging.error(f"Could not write a default config.toml: {exc}")
        if written is not None:
            logging.warning(f"No config.toml for this cruise: wrote {written} with DEFAULT SETTINGS.")
    config = ek80_ingestion.load_config_with_defaults(cruise_config_path, DEFAULTS_CONFIG_PATH)
    config_path = cruise_config_path if cruise_config_path.exists() else DEFAULTS_CONFIG_PATH
    unchecked = (config.get('CONFIG_CHECKED') is False) or not cruise_config_path.exists()
    if unchecked:
        config['_CONFIG_SOURCE'] = config['_CONFIG_SOURCE'] + ' - DEFAULT SETTINGS, NOT YET CHECKED'
        logging.warning("*** REVIEW REQUIRED: this build uses DEFAULT SETTINGS that have not been checked on this "
                        f"cruise's data. Check the seabed line on the echogram and the Calibration & QC tab in the dashboard, "
                        f"edit {cruise_config_path}, then set CONFIG_CHECKED = true. ***")
    if DEFAULTS_CONFIG_PATH.exists():
        config['_DEFAULTS_SHA256'] = file_sha256(DEFAULTS_CONFIG_PATH)
    config['_CONFIG_UNCHECKED'] = bool(unchecked)
    config['CRUISE_ID'] = cruise_id
    if args.range_bin is not None:
        config['MVBS_RANGE_BIN_M'] = args.range_bin
    if args.ping_bin is not None:
        config['MVBS_PING_TIME_BIN'] = args.ping_bin
    echogram_processing.validate_config(config)

    if not raw_dir.exists():
        logging.error(f"Raw directory not found: {raw_dir}")
        sys.exit(1)

    files = sorted(p for p in raw_dir.iterdir() if p.suffix.lower() == '.raw')
    if not files:
        logging.warning(f"No .raw files found in {raw_dir}")
        sys.exit(0)

    cruise_out = BASE_DIR / "processed" / cruise_id
    cruise_out.mkdir(parents=True, exist_ok=True)

    if not (cruise_dir / 'stations.csv').exists():
        write_stations_from_gps(cruise_dir, files)       # also writes nav.csv if the ship moved
    nav_path = cruise_dir / 'nav.csv'
    nav = echogram_analysis.load_nav(nav_path) if nav_path.exists() else None
    if nav is not None:
        logging.info(f"Read {len(nav)} navigation fixes from nav.csv")

    collected, file_rows, failures = {}, [], []
    for file_path in files:
        logging.info(f"Processing: {file_path.name}")
        try:
            tables, files_row = process_one_file(
                file_path, config, cruise_id, cruise_out, nav=nav, write_sv=not args.no_sv_store,
                write_nc=args.write_netcdf)
            for name, df in tables.items():
                collected.setdefault(name, []).append(df)
            file_rows.append(files_row)
        except Exception as exc:
            logging.error(f"Failed to process {file_path.name}: {exc}")
            logging.debug("Traceback for %s", file_path.name, exc_info=True)
            failures.append((file_path.name, str(exc)))
            continue

    if not file_rows:
        logging.warning("No data processed; database not updated.")
        for name, reason in failures:
            logging.warning(f"  {name}: {reason}")
        sys.exit(1)

    tables = {name: pd.concat(frames, ignore_index=True) for name, frames in collected.items()}
    for name in ('echogram_targets', 'echogram_schools', 'echogram_track'):
        if name in tables or name != 'echogram_track':      # track is optional; the others always exist (possibly empty)
            tables[name] = typed_table(name, tables.get(name, pd.DataFrame()))
    files_all = pd.DataFrame(file_rows)

    stations_path = cruise_dir / 'stations.csv'
    stations_df = load_stations(stations_path)
    if stations_df.empty:
        logging.warning(f"No stations.csv entries found at {stations_path}; positions will be empty.")
    files_all, stations_all = attach_positions(files_all, stations_df, cruise_id)
    tables['echogram_files'] = files_all
    tables['echogram_stations'] = stations_all

    # Optional Parquet copies beside the database, so other tools can read the tables directly.
    if args.write_parquet:
        for name, df in tables.items():
            df.to_parquet(cruise_out / f"{name.replace('echogram_', '')}.parquet", index=False)

    logging.info(f"Committing {len(tables['echogram_mvbs'])} MVBS cells, {len(tables['echogram_seabed'])} seabed pings, "
                 f"{len(tables['echogram_targets'])} targets and {len(tables['echogram_schools'])} school/layer "
                 f"detection(s) from {len(files_all)} file(s) to {args.db}")
    commit_to_duckdb(tables, cruise_id, args.db, config_path, config, build_started)

    if failures:
        logging.warning(f"Pipeline finished with {len(failures)} failed files:")
        for name, reason in failures:
            logging.warning(f"  {name}: {reason}")
    else:
        logging.info("Pipeline finished successfully.")
    if config.get('_CONFIG_UNCHECKED'):
        logging.warning("*** REVIEW REQUIRED: the settings used were DEFAULTS, not checked on this cruise. Open the dashboard, "
                        "look at the seabed line and the Calibration & QC tab, edit config.toml, then set CONFIG_CHECKED = true. ***")


if __name__ == "__main__":
    main()
