"""Echogram analysis helpers: block averaging, layer integration, frequency difference,
region statistics, school/layer detection and track interpolation.

All functions are plain NumPy/pandas so the build (main.py) and the dashboard
(echogram_holoviews.py) can share them and they can be tested on synthetic data.

Conventions: Sv is volume backscattering strength (dB re 1 m-1); sv is the same quantity
in linear units (m-1).  Averaging is always done in the linear domain.  sa is the area
backscattering coefficient (m2 m-2) = sum of sv x layer height; NASC (m2 nmi-2) = 4 pi 1852^2 x sa.
"""

import numpy as np
import pandas as pd
from scipy import ndimage

NASC_FACTOR = 4.0 * np.pi * 1852.0 ** 2          # converts sa (m2/m2) to NASC (m2/nmi2)
EARTH_RADIUS_M = 6371008.8


# ---------------------------------------------------------------------------------------------
# Block averaging (used for the smoothed display and for school detection)
# ---------------------------------------------------------------------------------------------
def block_average_sv(times, depth, sv_db, range_m, n_pings, min_fraction=0.5):
    """Averages Sv (dB) over blocks of range samples and pings, in the linear domain.

    ``sv_db`` has shape (range, ping).  A block counts only if at least ``min_fraction`` of its
    samples are valid, so mostly-removed areas stay blank.  Returns (times, depth, sv_db) of the
    block grid.  With range_m <= 0 and n_pings <= 1 the inputs are returned unchanged.
    """
    if range_m <= 0 and n_pings <= 1:
        return times, depth, sv_db
    dr = float(np.nanmedian(np.diff(depth)))
    n_r = max(1, int(round(range_m / dr))) if range_m > 0 else 1
    n_p = max(1, int(n_pings))
    nr, npg = sv_db.shape[0] // n_r * n_r, sv_db.shape[1] // n_p * n_p
    blocks = sv_db[:nr, :npg].reshape(nr // n_r, n_r, npg // n_p, n_p)
    valid = np.isfinite(blocks)
    with np.errstate(invalid='ignore', over='ignore'):
        power = np.nan_to_num(np.power(10.0, blocks / 10.0, dtype='float64'))
    count = valid.sum(axis=(1, 3))
    mean = power.sum(axis=(1, 3)) / np.maximum(count, 1)
    with np.errstate(divide='ignore'):
        out = 10 * np.log10(np.maximum(mean, 1e-30))
    out = np.where(count / (n_r * n_p) >= min_fraction, out, np.nan).astype('float32')
    new_depth = depth[:nr].reshape(-1, n_r).mean(axis=1)
    ticks = np.asarray(times)[:npg].astype('datetime64[ns]').astype('int64').reshape(-1, n_p).mean(axis=1)
    return ticks.astype('int64').astype('datetime64[ns]'), new_depth, out


# ---------------------------------------------------------------------------------------------
# Layers and integration (from the MVBS grid)
# ---------------------------------------------------------------------------------------------
def cell_height_m(depths):
    """Height of one MVBS cell, inferred from the depth values present."""
    d = np.unique(np.asarray(depths, dtype=float))
    return float(np.median(np.diff(d))) if len(d) > 1 else 1.0


def integrate_layers(mvbs, layer_m=10.0, time_bin=None):
    """Integrates MVBS cells into depth layers (and optionally coarser time bins).

    ``mvbs`` needs columns frequency_khz, ping_time, depth_m, mvbs_db (one row per non-empty cell).
    Returns one row per (frequency, time bin, layer) with the mean Sv, sa, NASC and the share of
    the layer's cells that held data (``coverage``).  Cells removed by cleaning are not counted.
    """
    cols = ['frequency_khz', 'time_bin', 'layer_top_m', 'layer_bottom_m', 'n_cells', 'coverage',
            'mean_sv_db', 'sa_m2_m2', 'nasc_m2_nmi2']
    df = mvbs.dropna(subset=['mvbs_db']).copy()
    if df.empty:
        return pd.DataFrame(columns=cols)
    dz = cell_height_m(df['depth_m'])
    df['layer_top_m'] = np.floor(df['depth_m'] / layer_m) * layer_m
    df['sv'] = np.power(10.0, df['mvbs_db'] / 10.0)
    df['time_bin'] = df['ping_time'].dt.floor(time_bin) if time_bin else df['ping_time']
    g = df.groupby(['frequency_khz', 'time_bin', 'layer_top_m'], sort=True)
    out = g.agg(n_cells=('sv', 'size'), sv_mean=('sv', 'mean'), sv_sum=('sv', 'sum')).reset_index()
    cells_per_layer = max(layer_m / dz, 1.0)
    n_time = (df.groupby('time_bin')['ping_time'].nunique().reindex(out['time_bin']).to_numpy()
              if time_bin else np.ones(len(out)))
    out['coverage'] = (out['n_cells'] / (cells_per_layer * n_time)).clip(upper=1.0)
    out['mean_sv_db'] = 10 * np.log10(out['sv_mean'])
    # sa per time step: sum of sv x height, averaged over the time steps inside the bin
    out['sa_m2_m2'] = out['sv_sum'] * dz / np.maximum(n_time, 1)
    out['nasc_m2_nmi2'] = out['sa_m2_m2'] * NASC_FACTOR
    out['layer_bottom_m'] = out['layer_top_m'] + layer_m
    return out[cols]


def frequency_difference(mvbs, f_high, f_low):
    """Sv(f_high) - Sv(f_low) per MVBS cell, where both frequencies have data (dB)."""
    a = mvbs[np.isclose(mvbs['frequency_khz'], f_high)].set_index(['ping_time', 'depth_m'])['mvbs_db']
    b = mvbs[np.isclose(mvbs['frequency_khz'], f_low)].set_index(['ping_time', 'depth_m'])['mvbs_db']
    diff = (a - b).dropna().rename('delta_db').reset_index()
    return diff


def region_stats(mvbs, t0, t1, d0, d1):
    """Statistics of the MVBS cells whose centres lie inside a time x depth box, per frequency.

    Also gives the dB difference between the highest and lowest frequency when there are two or more.
    """
    t0, t1 = pd.Timestamp(t0), pd.Timestamp(t1)
    d_lo, d_hi = sorted((float(d0), float(d1)))
    if mvbs.empty:
        return pd.DataFrame()
    dz = cell_height_m(mvbs['depth_m'])
    times = np.sort(mvbs['ping_time'].unique())
    dt = np.median(np.diff(times)) if len(times) > 1 else np.timedelta64(0, 'ns')
    centre_t = mvbs['ping_time'] + pd.Timedelta(dt) / 2
    centre_d = mvbs['depth_m'] + dz / 2
    inside = mvbs[(centre_t >= min(t0, t1)) & (centre_t <= max(t0, t1)) & (centre_d >= d_lo) & (centre_d <= d_hi)]
    rows = []
    for f, g in inside.groupby('frequency_khz'):
        sv = np.power(10.0, g['mvbs_db'].to_numpy() / 10.0)
        rows.append(dict(frequency_khz=f, n_cells=len(g), mean_sv_db=10 * np.log10(sv.mean()),
                         sa_m2_m2=float(sv.sum() * dz / max(g['ping_time'].nunique(), 1)),
                         nasc_m2_nmi2=float(sv.sum() * dz / max(g['ping_time'].nunique(), 1) * NASC_FACTOR)))
    out = pd.DataFrame(rows)
    if len(out) >= 2:
        out['delta_hi_minus_lo_db'] = out['mean_sv_db'].iloc[-1] - out['mean_sv_db'].iloc[0]
    return out


# ---------------------------------------------------------------------------------------------
# School / layer detection
# ---------------------------------------------------------------------------------------------
def detect_schools(times, depth, sv_db, threshold_db, min_thickness_m=2.0, min_length_s=5.0, max_gap_m=1.0,
                   max_gap_s=2.0):
    """Finds connected areas of backscatter above a threshold (schools or scattering layers).

    ``sv_db`` is a (range, ping) array, normally the block-averaged cleaned Sv.  Small gaps are
    closed, areas thinner than ``min_thickness_m`` or shorter than ``min_length_s`` are dropped.
    Returns one row per detection with its extent and backscatter.
    """
    cols = ['school_id', 'start_time', 'end_time', 'duration_s', 'depth_top_m', 'depth_bottom_m',
            'thickness_m', 'mean_depth_m', 'mean_sv_db', 'max_sv_db', 'n_cells']
    times = np.asarray(times).astype('datetime64[ns]')
    if sv_db.size == 0 or len(times) < 2:
        return pd.DataFrame(columns=cols)
    dz = float(np.nanmedian(np.diff(depth)))
    step_s = float(np.nanmedian(np.diff(times)) / np.timedelta64(1, 's'))
    mask = np.isfinite(sv_db) & (sv_db >= threshold_db)
    structure = np.ones((max(1, int(round(max_gap_m / dz)) * 2 + 1), max(1, int(round(max_gap_s / step_s)) * 2 + 1)), bool)
    closed = ndimage.binary_closing(mask, structure=structure) & (np.isfinite(sv_db) | mask)
    labels, n = ndimage.label(closed)
    rows = []
    for k, sl in enumerate(ndimage.find_objects(labels), start=1):
        if sl is None:
            continue
        sel = labels[sl] == k
        vals = sv_db[sl][sel & np.isfinite(sv_db[sl])]
        if vals.size == 0:
            continue
        r0, r1, c0, c1 = sl[0].start, sl[0].stop, sl[1].start, sl[1].stop
        thickness = (r1 - r0) * dz
        duration = (c1 - c0) * step_s
        if thickness < min_thickness_m or duration < min_length_s:
            continue
        rows.append(dict(
            school_id=len(rows) + 1, start_time=pd.Timestamp(times[c0]), end_time=pd.Timestamp(times[c1 - 1]),
            duration_s=duration, depth_top_m=float(depth[r0] - dz / 2), depth_bottom_m=float(depth[r1 - 1] + dz / 2),
            thickness_m=thickness, mean_depth_m=float(np.mean(depth[sl[0]][sel.any(axis=1)])),
            mean_sv_db=float(10 * np.log10(np.mean(np.power(10.0, vals / 10.0)))), max_sv_db=float(vals.max()),
            n_cells=int(sel.sum())))
    return pd.DataFrame(rows, columns=cols)


# ---------------------------------------------------------------------------------------------
# Navigation track
# ---------------------------------------------------------------------------------------------
def haversine_m(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi, dlmb = p2 - p1, np.radians(lon2) - np.radians(lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(a))


def load_nav(path):
    """Reads a navigation log: columns time_utc, lat, lon (decimal degrees).  Returns a frame sorted by time, or an empty frame."""
    cols = ['time_utc', 'lat', 'lon']
    try:
        df = pd.read_csv(path)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame(columns=cols)
    df = df.rename(columns={'time': 'time_utc', 'datetime': 'time_utc'})
    missing = [c for c in cols if c not in df]
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}; need time_utc, lat, lon")
    df['time_utc'] = pd.to_datetime(df['time_utc'], utc=True).dt.tz_convert(None)
    df[['lat', 'lon']] = df[['lat', 'lon']].apply(pd.to_numeric, errors='coerce')
    return df.dropna().sort_values('time_utc').reset_index(drop=True)[cols]


def interpolate_track(nav, ping_times, max_gap_s=120.0):
    """Interpolates a navigation log to ping times.  Returns lat, lon, distance_m (from the first
    ping) and speed_kn per ping; NaN where the log does not cover the ping within ``max_gap_s``."""
    pt = pd.to_datetime(np.asarray(ping_times))
    cols = ['ping_time', 'lat', 'lon', 'distance_m', 'speed_kn']
    if nav.empty:
        return pd.DataFrame(columns=cols)
    t_nav = nav['time_utc'].to_numpy().astype('datetime64[ns]').astype('int64') / 1e9
    t_png = pt.to_numpy().astype('datetime64[ns]').astype('int64') / 1e9
    lat = np.interp(t_png, t_nav, nav['lat'].to_numpy(), left=np.nan, right=np.nan)
    lon = np.interp(t_png, t_nav, nav['lon'].to_numpy(), left=np.nan, right=np.nan)
    nearest = np.abs(t_png[:, None] - t_nav[None, :]).min(axis=1) if len(t_nav) * len(t_png) < 5e7 else np.zeros(len(t_png))
    bad = nearest > max_gap_s
    lat[bad], lon[bad] = np.nan, np.nan
    step = haversine_m(lat[:-1], lon[:-1], lat[1:], lon[1:])
    dt = np.diff(t_png)
    dist = np.concatenate([[0.0], np.nancumsum(np.where(np.isfinite(step), step, 0.0))])
    speed = np.concatenate([[np.nan], step / np.maximum(dt, 1e-6) / 0.514444])
    speed = pd.Series(speed).rolling(5, min_periods=1, center=True).median().to_numpy()
    return pd.DataFrame({'ping_time': pt, 'lat': lat, 'lon': lon, 'distance_m': dist, 'speed_kn': speed})


# ---------------------------------------------------------------------------------------------
# Transient noise
# ---------------------------------------------------------------------------------------------
def transient_noise_mask(sv_db, depth, depth_bin_m=5.0, side_pings=25, exclude_above_m=10.0, threshold_db=12.0):
    """Flags short bursts of noise: depth bins whose level in one ping is far above the pings around it.

    Method after Ryan et al. (2015): ``sv_db`` is (range, ping).  Each ping's mean level (power mean) in
    every depth bin is compared with the median level of the ``side_pings`` pings before and after it;
    if it exceeds both by more than ``threshold_db`` the whole bin of that ping is flagged (at the
    start or end of a recording only the available side is used).  Bins shallower than
    ``exclude_above_m`` are never flagged.  Returns a boolean (range, ping) array, True = transient noise.
    """
    n_range, n_ping = sv_db.shape
    dz = float(np.nanmedian(np.diff(depth)))
    n_bin = max(1, int(round(depth_bin_m / dz)))
    n_blocks = n_range // n_bin
    flagged = np.zeros(sv_db.shape, bool)
    if n_blocks == 0 or n_ping < 3:
        return flagged
    blocks = sv_db[:n_blocks * n_bin].reshape(n_blocks, n_bin, n_ping)
    with np.errstate(invalid='ignore', over='ignore'):
        power = np.where(np.isfinite(blocks), np.power(10.0, blocks / 10.0), np.nan)
        level = 10 * np.log10(np.nanmean(power, axis=1))                 # (block, ping)
    left = np.full_like(level, np.nan)
    right = np.full_like(level, np.nan)
    for p in range(n_ping):
        lo, hi = max(0, p - side_pings), min(n_ping, p + side_pings + 1)
        if p > lo:
            left[:, p] = np.nanmedian(level[:, lo:p], axis=1)
        if hi > p + 1:
            right[:, p] = np.nanmedian(level[:, p + 1:hi], axis=1)
    with np.errstate(invalid='ignore'):
        over_left = (level - left) > threshold_db
        over_right = (level - right) > threshold_db
        has_left, has_right = np.isfinite(left), np.isfinite(right)
        hit = np.where(has_left & has_right, over_left & over_right, np.where(has_left, over_left, over_right))
    bin_depth = depth[:n_blocks * n_bin].reshape(n_blocks, n_bin).mean(axis=1)
    hit &= (bin_depth >= exclude_above_m)[:, None]
    hit &= np.isfinite(level)
    flagged[:n_blocks * n_bin] = np.repeat(hit, n_bin, axis=0)
    return flagged
