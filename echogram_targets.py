"""Single-target detection on target strength (TS) data.

Looks for isolated echoes from individual animals along each ping, using the usual criteria:
a peak in TS above a threshold, enough signal above the background noise, an echo length close to
the transmitted pulse length, and (for split-beam channels) a position in the beam close enough
to the axis that the beam-pattern correction is small.

The detector is a plain NumPy function so it can be tested on synthetic echoes.
"""

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

TARGET_COLUMNS = ['frequency_khz', 'ping_index', 'ping_time', 'depth_m', 'ts_db', 'ts_comp_db', 'snr_db', 'peak_over_background_db',
                  'along_deg', 'athwart_deg', 'beam_comp_db', 'echo_length_ratio', 'compensated']


def beam_compensation_db(along_deg, athwart_deg, bw_along_deg, bw_athwart_deg):
    """Two-way beam-pattern loss (dB, >= 0) for a target at the given split-beam angles.

    Uses the standard parabolic approximation: 6.0206 x [(2a/bw_a)^2 + (2b/bw_b)^2], where bw is the
    half-power (-3 dB) full beam width.  A target at the half-power angle gives 6 dB.
    """
    return 6.0206 * ((2.0 * along_deg / bw_along_deg) ** 2 + (2.0 * athwart_deg / bw_athwart_deg) ** 2)


def _cumulative_power(ts):
    """Running sums of linear power and of valid-sample counts along range, with a leading zero column."""
    valid = np.isfinite(ts)
    power = np.where(valid, np.power(10.0, np.where(valid, ts, 0.0) / 10.0), 0.0)
    zeros = np.zeros((ts.shape[0], 1))
    return (np.concatenate([zeros, np.cumsum(power, axis=1)], axis=1),
            np.concatenate([zeros, np.cumsum(valid.astype(float), axis=1)], axis=1))


def background_db(cum_power, cum_count, p, k, half_window, exclude):
    """Average level (dB, power mean) of the samples within ``half_window`` of sample k in ping p,
    leaving out ``exclude`` samples either side of k (the target's own echo).  NaN if too few samples."""
    n = cum_power.shape[1] - 1
    a, b = max(0, k - half_window), min(n - 1, k + half_window)
    ea, eb = max(0, k - exclude), min(n - 1, k + exclude)
    total = cum_power[p, b + 1] - cum_power[p, a] - (cum_power[p, eb + 1] - cum_power[p, ea])
    count = cum_count[p, b + 1] - cum_count[p, a] - (cum_count[p, eb + 1] - cum_count[p, ea])
    if count < 5 or total <= 0:
        return np.nan
    return 10.0 * np.log10(total / count)


def detect_single_targets(ts, sv, sv_noise, depth, ping_times, pulse_samples, frequency_khz,
                          along=None, athwart=None, bw_along=None, bw_athwart=None,
                          ts_min_db=-60.0, snr_min_db=10.0, min_length=0.8, max_length=1.8,
                          max_beam_comp_db=6.0, min_depth_m=5.0, max_depth_m=None,
                          background_window_m=5.0, min_peak_over_background_db=15.0):
    """Detects single targets in one channel.

    ts, sv, sv_noise: (ping, range) arrays in dB (target strength, calibrated Sv, background-noise Sv).
    depth: (range,) depth of each sample in metres.  pulse_samples: pulse length in samples.
    along, athwart: optional (ping, range) split-beam angles in degrees, with bw_along / bw_athwart the
    beam widths.  ``max_depth_m`` is a scalar or a (ping,) array (e.g. the seabed minus a margin).
    A peak must also stand at least ``min_peak_over_background_db`` above the average TS of the samples
    within ``background_window_m`` of it in the same ping, leaving out one pulse length either side of the
    peak (so the target does not raise its own background).  This is what keeps noise from being counted.
    Returns a DataFrame of detections (columns in TARGET_COLUMNS).
    """
    n_ping, n_range = ts.shape
    have_angles = along is not None and athwart is not None and bw_along and bw_athwart
    max_depth = np.full(n_ping, np.inf) if max_depth_m is None else np.broadcast_to(np.asarray(max_depth_m, float), (n_ping,))
    distance = max(1, int(round(pulse_samples)))
    dz = float(np.nanmedian(np.diff(depth[np.isfinite(depth)]))) if np.isfinite(depth).sum() > 1 else 1.0
    cum_power, cum_count = _cumulative_power(ts)
    half_window = max(3, int(round(background_window_m / dz / 2)))
    rows = []
    for p in range(n_ping):
        row = ts[p]
        usable = np.isfinite(row) & (depth >= min_depth_m) & (depth <= max_depth[p])
        if not usable.any():
            continue
        work = np.where(usable, row, -np.inf)
        peaks, _ = find_peaks(work, height=ts_min_db, distance=distance)
        for k in peaks:
            snr = sv[p, k] - sv_noise[p, k]
            if not np.isfinite(snr) or snr < snr_min_db:
                continue
            bg = background_db(cum_power, cum_count, p, k, half_window, distance)
            if not row[k] - bg >= min_peak_over_background_db:
                continue
            level = row[k] - 6.0                                  # echo length measured 6 dB below the peak
            lo = k
            while lo > 0 and np.isfinite(row[lo - 1]) and row[lo - 1] >= level:
                lo -= 1
            hi = k
            while hi < n_range - 1 and np.isfinite(row[hi + 1]) and row[hi + 1] >= level:
                hi += 1
            ratio = (hi - lo + 1) / pulse_samples
            if not (min_length <= ratio <= max_length):
                continue
            a = b = comp = np.nan
            if have_angles:
                w = np.power(10.0, row[lo:hi + 1] / 10.0)          # power-weighted angle over the echo
                a = float(np.average(along[p, lo:hi + 1], weights=w))
                b = float(np.average(athwart[p, lo:hi + 1], weights=w))
                comp = float(beam_compensation_db(a, b, bw_along, bw_athwart))
                if not np.isfinite(comp) or comp > max_beam_comp_db:
                    continue
            rows.append(dict(
                frequency_khz=float(frequency_khz), ping_index=p, ping_time=pd.Timestamp(ping_times[p]),
                depth_m=float(depth[k]), ts_db=float(row[k]), ts_comp_db=float(row[k] + (comp if have_angles else 0.0)),
                snr_db=float(snr), peak_over_background_db=float(row[k] - bg), along_deg=a, athwart_deg=b, beam_comp_db=comp,
                echo_length_ratio=float(ratio), compensated=bool(have_angles)))
    return pd.DataFrame(rows, columns=TARGET_COLUMNS)
