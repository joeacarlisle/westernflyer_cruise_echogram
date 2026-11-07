"""Western Flyer EK80 echogram dashboard (Panel + HoloViews + Datashader).

Run with:  python -m panel serve echogram_holoviews.py --show

Reads what main.py wrote:
    processed/wf_echogram.duckdb                 MVBS grid, seabed line, file list, build provenance
    processed/<cruise>/sv/<file_key>.zarr        cleaned full-resolution Sv (read lazily per file)

Large data is never sent to the browser cell by cell: Datashader (via
HoloViews ``rasterize``) turns the Sv array into an image sized to the screen
and re-renders it as you zoom.
"""

import functools
import io
import json
import pathlib
import sys

import duckdb
import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn
import xarray as xr
from bokeh.models import CustomJS
from holoviews.element.tiles import Tiles
from holoviews.operation.datashader import rasterize

BASE_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(BASE_DIR))          # so the dashboard finds echogram_analysis however it is launched
import echogram_analysis as analysis       # noqa: E402

# 1. ENGINE INITIALIZATION
pn.extension('tabulator')
hv.extension('bokeh')

DB_PATH = BASE_DIR / 'processed' / 'wf_echogram.duckdb'


# EXPORT STYLING (same approach as the CTD dashboard): the toolbar's Save button exports a PNG, and
# the dashboard theme makes the figure background transparent, which shows as a checkerboard.
def _white_figure(fig, background=True):
    """Opaque white figure, so exported PNGs have no transparent area."""
    if background:
        fig.background_fill_color = 'white'
        fig.background_fill_alpha = 0.999
    fig.border_fill_color = 'white'
    # The theme sets border_fill_alpha = 0.  Bokeh does not send a property equal to its default (1.0),
    # so the theme would win; 0.999 is sent and still exports fully opaque.
    fig.border_fill_alpha = 0.999


def _white_background(plot, element):
    """HoloViews hook: opaque white figure with dark labels and title."""
    fig = plot.state
    _white_figure(fig)
    for ax in list(fig.xaxis) + list(fig.yaxis):
        ax.axis_label_text_color = '#222222'
        ax.major_label_text_color = '#222222'
    if fig.title is not None:
        fig.title.text_color = '#222222'


def _white_border(plot, element):
    """Like _white_background but keeps the plot area's own colour (used for the map)."""
    _white_figure(plot.state, background=False)


# Stacked frequency panels export as one grid.  Its margins and toolbar strip are not part of any single
# figure, and at display scaling above 100 % part of the backing is unpainted, so the export has
# transparent areas.  This one-time patch paints white behind everything in a grid export.
_WHITE_EXPORT_JS = """
const M = Bokeh.Models.get('GridPlot');
const V = M && M.prototype.default_view;
if (V && !V.prototype.__white_export) {
  V.prototype.__white_export = true;
  const orig = V.prototype.export;
  V.prototype.export = function (...args) {
    const out = orig.apply(this, args);
    try {
      const c = out.canvas, ctx = out.ctx;
      if (c instanceof HTMLCanvasElement) {
        ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0);
        ctx.globalCompositeOperation = 'destination-over';
        ctx.fillStyle = 'white'; ctx.fillRect(0, 0, c.width, c.height); ctx.restore();
      }
    } catch (e) {}
    return out;
  };
}
"""
if pn.state.curdoc is not None:
    pn.state.curdoc.js_on_event('document_ready', CustomJS(code=_WHITE_EXPORT_JS))

VIEW_FULL = 'Full resolution (cleaned Sv)'
VIEW_MVBS = 'MVBS grid (gridded mean)'
VIEW_DIFF = 'Frequency difference (MVBS)'
VIEW_RAW = 'Uncalibrated echo level'


# 2. DATA ACCESS
def query_df(sql, params=None):
    """Runs a read-only query against the echogram DuckDB; returns an empty frame if it is missing."""
    if not DB_PATH.exists():
        return pd.DataFrame()
    with duckdb.connect(str(DB_PATH), read_only=True) as con:
        return con.execute(sql, params or []).df()


@functools.lru_cache(maxsize=6)
def load_sv_channel(zarr_rel_path, freq_khz):
    """Loads one channel of the cleaned Sv store: (ping times, depth axis m, Sv[range, ping]).

    Cached so switching colour limits or the seabed overlay does not re-read the file.
    """
    ds = xr.open_zarr(BASE_DIR / zarr_rel_path, consolidated=False)
    idx = int(np.argmin(np.abs(ds['frequency_khz'].values - freq_khz)))
    depth = ds['display_depth_m'].isel(channel=idx).values
    valid = np.isfinite(depth)
    sv = ds['Sv'].isel(channel=idx).values[:, valid]          # (ping, range)
    return ds['ping_time'].values, depth[valid], sv.T.astype('float32')


SMOOTHING = {                       # label -> (range block in m, pings per block)
    'Off': (0.0, 1),
    'Light (0.25 m x 3 pings)': (0.25, 3),
    'Medium (0.5 m x 5 pings)': (0.5, 5),
    'Strong (1 m x 9 pings)': (1.0, 9),
}
SMOOTH_MIN_FRACTION = 0.5           # a block with less than this share of valid samples stays blank


@functools.lru_cache(maxsize=6)
def load_sv_smoothed(zarr_rel_path, freq_khz, range_m, n_pings):
    """Cleaned Sv averaged over blocks of samples and pings, for a smoother picture.

    Averaging is done on linear power (not on dB).  A block counts only if at least half of its
    samples survived cleaning, so mostly-removed areas stay blank instead of being filled in
    from a few leftover strong samples.
    """
    times, depth, sv = load_sv_channel(zarr_rel_path, freq_khz)
    return analysis.block_average_sv(times, depth, sv, range_m, n_pings, SMOOTH_MIN_FRACTION)


@functools.lru_cache(maxsize=2)
def load_raw_level(zarr_rel_path):
    """Uncalibrated echo level for every channel: (ping times, [(depth axis m, level_db[range, ping]), ...])."""
    ds = xr.open_zarr(BASE_DIR / zarr_rel_path, consolidated=False)
    panels = []
    for i in range(ds.sizes['channel']):
        depth = ds['depth_m'].isel(channel=i).values
        valid = np.isfinite(depth)
        level = ds['level_db'].isel(channel=i).values[:, valid]          # (ping, range)
        panels.append((depth[valid], level.T.astype('float32')))
    return ds['ping_time'].values, panels


def query_safe(sql, params=None):
    """Like query_df, but an empty frame when the table does not exist yet (older database)."""
    try:
        return query_df(sql, params)
    except duckdb.CatalogException:
        return pd.DataFrame()


# 3. WIDGETS
files_all = query_df("SELECT * FROM echogram_files ORDER BY cruise_id, first_ping") \
    if DB_PATH.exists() else pd.DataFrame()
cruises = sorted(files_all['cruise_id'].unique().tolist()) if not files_all.empty else []

cruise_select = pn.widgets.Select(name='Cruise', options=cruises, value=cruises[0] if cruises else None)
file_select = pn.widgets.Select(name='Recording (.raw file)')
view_select = pn.widgets.RadioButtonGroup(name='View', options=[VIEW_FULL, VIEW_MVBS, VIEW_DIFF, VIEW_RAW], value=VIEW_FULL,
                                          button_type='default', orientation='vertical')
depth_max = pn.widgets.IntSlider(name='Maximum depth shown (m)', start=10, end=1000, value=130, step=5)
plot_height = pn.widgets.IntSlider(name='Plot height per frequency (px)', start=250, end=1000, value=380, step=10)
seabed_toggle = pn.widgets.Checkbox(name='Show detected seabed', value=True)
smooth_select = pn.widgets.Select(name='Smoothing (full resolution view)', options=list(SMOOTHING), value='Off')
lock_clim = pn.widgets.Checkbox(name='Lock colour ranges (compare recordings)', value=False)
sv_threshold = pn.widgets.IntSlider(name='Hide Sv below (dB; -110 = off)', start=-110, end=-30, value=-110, step=1)
diff_clim = pn.widgets.RangeSlider(name='Difference range (dB; difference view)', start=-80, end=80, value=(-20, 20), step=1)
layers_toggle = pn.widgets.Checkbox(name='Show layer boundaries', value=False)
layer_select = pn.widgets.Select(name='Layer thickness (m)', options=[5, 10, 20, 50], value=10)
targets_toggle = pn.widgets.Checkbox(name='Show single targets', value=True)
schools_toggle = pn.widgets.Checkbox(name='Show detected schools / layers', value=True)

# One colour-range slider per frequency found in the selected recording.  They are rebuilt whenever the
# recording changes, so the dashboard follows whatever frequencies the file contains.
freq_sliders = {}                 # frequency (kHz) -> RangeSlider
sliders_box = pn.Column()         # sidebar container for those sliders
clim_version = pn.widgets.IntInput(value=0)    # bumped when any colour range changes, to redraw (not shown)
locked_ranges = {}                # frequency (kHz) -> (low, high) dB, remembered while the lock is on
_setting_slider = False           # True while the code, not the user, moves a colour slider


def _bump():
    clim_version.value += 1


def _file_row(file_key):
    sub = files_all[files_all['file_key'] == file_key]
    return sub.iloc[0] if len(sub) else None


def _file_freqs(row):
    """Frequencies (kHz) detected in a recording, in channel order."""
    return [float(f) for f in str(row['frequencies_khz']).split(',') if f]


def _make_slider(freq):
    slider = pn.widgets.RangeSlider(name=f"Colour range {freq:g} kHz (dB)",
                                    start=-110, end=0, value=(-85, -40), step=1)

    def on_change(*events):
        # A drag can deliver 'value' and 'value_throttled' together, so accept any number of events.
        if _setting_slider or not events:
            return
        if lock_clim.value:
            locked_ranges[freq] = tuple(events[-1].new)      # manual change while locked is remembered
        if any(e.name == 'value' for e in events):
            _bump()

    slider.param.watch(on_change, ['value', 'value_throttled'])    # drags and key presses
    return slider


def _set_slider(slider, value):
    """Moves a colour slider from code without it counting as a manual change."""
    global _setting_slider
    _setting_slider = True
    try:
        slider.value = tuple(value)
    finally:
        _setting_slider = False


def _update_files(*_):
    sub = files_all[files_all['cruise_id'] == cruise_select.value]
    file_select.options = {r.filename: r.file_key for r in sub.itertuples()}
    if len(sub):
        file_select.value = sub['file_key'].iloc[0]


def _update_freqs(*_):
    """Rebuilds the colour sliders for the frequencies in the selected recording."""
    row = _file_row(file_select.value)
    freq_sliders.clear()
    if row is None:
        sliders_box[:] = []
        return
    freq_sliders.update({f: _make_slider(f) for f in _file_freqs(row)})
    sliders_box[:] = list(freq_sliders.values())
    if np.isfinite(row['seabed_max_m']):
        depth_max.value = int(max(10, np.ceil(row['seabed_max_m'] * 1.4 / 5) * 5))


def _auto_range(row, freq, view):
    """Colour range for one frequency from the data (1st to 99.5th percentile), or None if there is no data."""
    try:
        if view == VIEW_FULL and row['sv_zarr_path']:
            _, _, sv = load_sv_smoothed(row['sv_zarr_path'], float(freq), *SMOOTHING[smooth_select.value])
            vals = sv[::max(1, sv.shape[0] // 600), ::2]
        else:
            vals = query_df("SELECT mvbs_db FROM echogram_mvbs WHERE file_key = ? AND frequency_khz = ?",
                            [row['file_key'], float(freq)])['mvbs_db'].to_numpy()
        vals = vals[np.isfinite(vals)]
        if vals.size < 10:
            return None
        lo, hi = np.percentile(vals, [1, 99.5])
    except Exception:
        return None
    lo, hi = float(np.floor(lo)), float(np.ceil(hi))
    if hi - lo < 10:
        hi = lo + 10
    return max(lo, -110.0), min(hi, 0.0)


def _auto_diff_range(row):
    """Colour range for the difference view from the data (1st to 99th percentile), or None."""
    try:
        freqs = _file_freqs(row)
        if len(freqs) < 2:
            return None
        mv = query_df("SELECT frequency_khz, ping_time, depth_m, mvbs_db FROM echogram_mvbs WHERE file_key = ?", [row['file_key']])
        d = analysis.frequency_difference(mv, max(freqs), min(freqs))['delta_db'].to_numpy()
        d = d[np.isfinite(d)]
        if d.size < 10:
            return None
        lo, hi = np.percentile(d, [1, 99])
    except Exception:
        return None
    lo, hi = int(np.floor(lo)), int(np.ceil(hi))
    if hi - lo < 10:
        lo, hi = lo - 5, hi + 5
    return max(lo, -80), min(hi, 80)


def _auto_clim(*_):
    """Sets every frequency's colour range from the data, or from the lock if it is on.

    Different frequencies have very different Sv levels, so each one gets its own range.
    The sliders can still be moved by hand afterwards.
    """
    row = _file_row(file_select.value)
    view = view_select.value
    smooth_select.disabled = (view != VIEW_FULL)    # smoothing only applies to the full-resolution view
    sv_threshold.disabled = view in (VIEW_RAW, VIEW_DIFF)
    diff_clim.disabled = (view != VIEW_DIFF)
    for slider in freq_sliders.values():
        slider.disabled = view in (VIEW_RAW, VIEW_DIFF)   # those views set their own colour scale
    if row is not None and view == VIEW_DIFF and not lock_clim.value:
        rng = _auto_diff_range(row)                      # fit the difference scale to this recording
        if rng is not None:
            diff_clim.value = rng
    if row is None or view in (VIEW_RAW, VIEW_DIFF):
        _bump()
        return
    for freq, slider in freq_sliders.items():
        if lock_clim.value and freq in locked_ranges:    # locked: reuse this frequency's remembered range
            _set_slider(slider, locked_ranges[freq])
            continue
        rng = _auto_range(row, freq, view)
        if rng is None:
            continue
        _set_slider(slider, rng)
        if lock_clim.value:
            locked_ranges[freq] = rng                    # first time this frequency is seen while locked
    _bump()


def _lock_toggled(event):
    if event.new:
        locked_ranges.update({f: tuple(sl.value) for f, sl in freq_sliders.items()})
    else:
        locked_ranges.clear()
        _auto_clim()                                     # unlocked: go back to automatic ranges


def _reset_display(*_):
    """Puts the display controls back to their starting values.

    Colour ranges return to automatic, the lock is released and the thresholds, smoothing, depth
    and plot height go back to their defaults.  The cruise, recording, view and the overlay
    checkboxes are left as they are.
    """
    locked_ranges.clear()
    lock_clim.value = False                  # also triggers the automatic colour ranges
    sv_threshold.value = -110
    smooth_select.value = 'Off'
    plot_height.value = 380
    row = _file_row(file_select.value)
    depth_max.value = int(max(10, np.ceil(row['seabed_max_m'] * 1.4 / 5) * 5)) \
        if row is not None and np.isfinite(row['seabed_max_m']) else 130
    _auto_clim()                             # also refits the difference range when that view is shown
    if view_select.value != VIEW_DIFF:
        diff_clim.value = (-20, 20)


reset_button = pn.widgets.Button(name='Reset display', button_type='warning', width=290)
reset_button.on_click(_reset_display)

cruise_select.param.watch(_update_files, 'value')
file_select.param.watch(_update_freqs, 'value')
file_select.param.watch(_auto_clim, 'value')      # registered after _update_freqs, so it runs second
view_select.param.watch(_auto_clim, 'value')
smooth_select.param.watch(_auto_clim, 'value')    # smoothing narrows the spread of values, so re-fit the ranges
lock_clim.param.watch(_lock_toggled, 'value')
_update_files()
_update_freqs()
_auto_clim()


# 4. ECHOGRAM VIEW
def _seabed_curve(file_key):
    sb = query_df("SELECT ping_time, seabed_depth_m FROM echogram_seabed WHERE file_key = ? ORDER BY ping_time",
                  [file_key])
    return hv.Curve(sb, kdims=['ping_time'], vdims=['seabed_depth_m']).opts(
        color='red', line_width=1.5, tools=[])


def _targets_overlay(file_key, freq):
    t = query_df("SELECT ping_time, depth_m, ts_comp_db, snr_db, peak_over_background_db FROM echogram_targets "
                 "WHERE file_key = ? AND frequency_khz = ?", [file_key, float(freq)])
    if t.empty:
        return None
    return hv.Points(t, kdims=['ping_time', 'depth_m'], vdims=['ts_comp_db', 'snr_db', 'peak_over_background_db']).opts(
        marker='circle', size=12, fill_alpha=0, line_color='red', line_width=2, tools=['hover'])


def _schools_overlay(file_key, freq):
    s = query_df("SELECT school_id, start_time, end_time, depth_top_m, depth_bottom_m, mean_sv_db, thickness_m "
                 "FROM echogram_schools WHERE file_key = ? AND frequency_khz = ?", [file_key, float(freq)])
    if s.empty:
        return None
    data = [(r.start_time, r.depth_top_m, r.end_time, r.depth_bottom_m, int(r.school_id), float(r.mean_sv_db), float(r.thickness_m))
            for r in s.itertuples()]
    return hv.Rectangles(data, kdims=['x0', 'y0', 'x1', 'y1'], vdims=['school', 'mean_sv_db', 'thickness_m']).opts(
        fill_alpha=0, line_color='white', line_dash='dashed', line_width=2, tools=['hover'])


def _layer_lines(ymax, layer_m):
    lines = [hv.HLine(float(y)).opts(color='white', line_dash='dotted', line_width=1, alpha=0.7)
             for y in np.arange(layer_m, ymax, layer_m)]
    return hv.Overlay(lines) if lines else None


def _time_range_hook(start_ms, end_ms, plot, element):
    """Sets the time axis to the recording, always earliest to latest (never reversed or stretched).

    Bokeh stores datetime axes as milliseconds since 1970, so this does not depend on how a
    particular HoloViews version converts a datetime axis limit.
    """
    rng = plot.state.x_range
    rng.start, rng.end = float(start_ms), float(end_ms)


def _decorate(plot, file_key, freq, ymax, show_seabed, show_layers, layer_m, show_targets, show_schools):
    """Adds the seabed line, layer boundaries, single targets and detected schools to one panel."""
    hv.streams.Tap(source=plot, x=None, y=None, subscribers=[_on_tap])
    if show_seabed:
        plot = plot * _seabed_curve(file_key)
    if show_layers and (lines := _layer_lines(ymax, float(layer_m))) is not None:
        plot = plot * lines
    if show_schools and (rects := _schools_overlay(file_key, freq)) is not None:
        plot = plot * rects
    if show_targets and (pts := _targets_overlay(file_key, freq)) is not None:
        plot = plot * pts
    out = plot * hv.DynamicMap(_cursor_lines, streams=[cursor])
    hooks = [_white_background]
    row = _file_row(file_key)
    if row is not None and pd.notna(row['first_ping']) and pd.notna(row['last_ping']):
        t_a, t_b = sorted((pd.Timestamp(row['first_ping']), pd.Timestamp(row['last_ping'])))
        hooks.append(functools.partial(_time_range_hook, t_a.value // 1_000_000, t_b.value // 1_000_000))
    return out.opts(hooks=hooks)


Cursor = hv.streams.Stream.define('Cursor', x=None, y=None)
cursor = Cursor()                     # the clicked position, shared by every panel
details_pane = pn.pane.Markdown('Click anywhere on an echogram to read the values at that point '
                                '(a cyan cross-hair appears on every frequency).', sizing_mode='stretch_width')


def _to_datetime64(x):
    """A click on a time axis may arrive as datetime64, a pandas Timestamp or milliseconds since 1970."""
    if isinstance(x, (int, float, np.integer, np.floating)):
        return np.datetime64(int(x), 'ms').astype('datetime64[ns]')
    return np.datetime64(pd.Timestamp(x)).astype('datetime64[ns]')


def _cursor_lines(x, y):
    """Cross-hair at the clicked position (invisible until something has been clicked)."""
    shown = x is not None
    when = _to_datetime64(x) if shown else np.datetime64('2000-01-01T00:00:00', 'ns')
    v = hv.VLine(when).opts(color='cyan', line_width=1.5, alpha=0.9 if shown else 0.0)
    h = hv.HLine(float(y) if shown else 0.0).opts(color='cyan', line_width=1, alpha=0.6 if shown else 0.0)
    return v * h


def _nearest(values, target):
    i = int(np.searchsorted(values, target))
    if i > 0 and (i == len(values) or abs(values[i - 1] - target) <= abs(values[i] - target)):
        i -= 1
    return min(i, len(values) - 1)


def _details_markdown(x, y):
    row = _file_row(file_select.value)
    if row is None:
        return ''
    t, depth_clicked = _to_datetime64(x), float(y)
    out = [f"**{pd.Timestamp(t):%Y-%m-%d %H:%M:%S} UTC**, depth **{depth_clicked:.1f} m**"]
    lines = ['| Frequency | Ping # | Cleaned Sv (dB) | Uncalibrated level (dB) | Seabed (m) |', '|---|---|---|---|---|']
    freqs = _file_freqs(row)
    raw = load_raw_level(row['raw_level_zarr_path']) if row['raw_level_zarr_path'] else None
    for k, f in enumerate(freqs):
        sv_txt = raw_txt = '-'
        ping = '-'
        if row['sv_zarr_path']:
            times, depth, sv = load_sv_channel(row['sv_zarr_path'], float(f))
            if len(times) and times[0] <= t <= times[-1] + np.timedelta64(5, 's'):
                ip, iz = _nearest(times, t), _nearest(depth, depth_clicked)
                ping = str(ip)
                sv_txt = f"{sv[iz, ip]:.1f}" if np.isfinite(sv[iz, ip]) else 'removed by cleaning'
        if raw is not None and k < len(raw[1]):
            rt, panels = raw
            rdepth, rlevel = panels[k]
            ip, iz = _nearest(rt, t), _nearest(rdepth, depth_clicked)
            raw_txt = f"{rlevel[iz, ip]:.1f}" if np.isfinite(rlevel[iz, ip]) else '-'
        sb = query_safe("SELECT seabed_depth_m FROM echogram_seabed WHERE file_key = ? AND frequency_khz = ? "
                        "ORDER BY abs(epoch(ping_time) - epoch(?::TIMESTAMP)) LIMIT 1", [row['file_key'], float(f), str(pd.Timestamp(t))])
        sb_txt = f"{sb['seabed_depth_m'].iloc[0]:.1f}" if not sb.empty and np.isfinite(sb['seabed_depth_m'].iloc[0]) else '-'
        lines.append(f"| {f:g} kHz | {ping} | {sv_txt} | {raw_txt} | {sb_txt} |")
    out.append('\n'.join(lines))
    tstr = str(pd.Timestamp(t))
    near = query_safe("SELECT frequency_khz, depth_m, ts_comp_db FROM echogram_targets WHERE file_key = ? "
                      "AND abs(epoch(ping_time) - epoch(?::TIMESTAMP)) < 5 AND abs(depth_m - ?) < 3", [row['file_key'], tstr, depth_clicked])
    if not near.empty:
        out.append("Single target nearby: " + "; ".join(f"{r.frequency_khz:g} kHz at {r.depth_m:.1f} m, TS {r.ts_comp_db:.1f} dB" for r in near.itertuples()))
    sch = query_safe("SELECT frequency_khz, school_id, mean_sv_db FROM echogram_schools WHERE file_key = ? AND start_time <= ?::TIMESTAMP "
                     "AND end_time >= ?::TIMESTAMP AND depth_top_m <= ? AND depth_bottom_m >= ?",
                     [row['file_key'], tstr, tstr, depth_clicked, depth_clicked])
    if not sch.empty:
        out.append("Inside detected school/layer: " + "; ".join(f"{r.frequency_khz:g} kHz #{r.school_id} (mean Sv {r.mean_sv_db:.1f} dB)" for r in sch.itertuples()))
    trk = query_safe("SELECT lat, lon, speed_kn, distance_m FROM echogram_track WHERE file_key = ? "
                     "ORDER BY abs(epoch(ping_time) - epoch(?::TIMESTAMP)) LIMIT 1", [row['file_key'], tstr])
    if not trk.empty and np.isfinite(trk['lat'].iloc[0]):
        r = trk.iloc[0]
        out.append(f"Ship: {r.lat:.4f}, {r.lon:.4f}, {r.speed_kn:.1f} kn, {r.distance_m / 1852:.2f} nmi from start of file")
    return '  \n'.join(out)


def _on_tap(x, y):
    if x is None or y is None:
        return
    try:
        details_pane.object = _details_markdown(x, y)
    except Exception as exc:          # a click must never break the dashboard
        details_pane.object = f"Could not read values at that point: {exc}"
    cursor.event(x=x, y=y)


def _stack(plots):
    """All frequencies in one column that shares its axes; scrolls if taller than the window."""
    layout = hv.Layout(plots).cols(1).opts(shared_axes=True)
    return pn.pane.HoloViews(layout, sizing_mode='stretch_width')


def _raw_level_impl(row, ymax, height, deco):
    """Every frequency stacked, as recorded: no calibration, noise removal or seabed masking."""
    if not row['raw_level_zarr_path']:
        return pn.pane.Alert('No uncalibrated echo level for this file (built with --no-sv-store, or the '
                             'preview was skipped). Rebuild to create it.', alert_type='warning')
    times, panels = load_raw_level(row['raw_level_zarr_path'])
    freqs = _file_freqs(row)
    plots = []
    for k, ((depth, level), freq) in enumerate(zip(panels, freqs)):
        # Each panel is scaled on its own: 1st to 99.8th percentile below the near-field strip.
        body = level[depth > 5]
        vmin, vmax = np.nanpercentile(body, [1, 99.8])
        quad = hv.QuadMesh((times, depth, level), kdims=['ping_time', 'depth_m'], vdims=['level'])
        img = rasterize(quad, aggregator='mean', dynamic=True).opts(
            cmap='viridis', clim=(float(vmin), float(vmax)), colorbar=True, invert_yaxis=True,
            ylim=(0, ymax), responsive=True, height=height, active_tools=['pan'], title=f"{freq:g} kHz",
            tools=['hover', 'wheel_zoom'],
            xlabel='Ping time (UTC)' if k == len(panels) - 1 else '', ylabel='Depth (m)',
            clabel='Echo level (dB, uncalibrated)')
        plots.append(_decorate(img, row['file_key'], freq, ymax, **deco))
    return _stack(plots)


def _calibrated_impl(row, view, ymax, height, smooth, threshold, deco):
    """Every frequency stacked, calibrated Sv (full resolution or MVBS grid)."""
    freqs = _file_freqs(row)
    if view == VIEW_FULL and not row['sv_zarr_path']:
        return pn.pane.Alert('This file was built with --no-sv-store; choose the MVBS view.',
                             alert_type='warning')
    mv_all = query_df("SELECT frequency_khz, ping_time, depth_m, mvbs_db FROM echogram_mvbs WHERE file_key = ?",
                      [row['file_key']]) if view == VIEW_MVBS else None
    plots = []
    for k, freq in enumerate(freqs):
        if view == VIEW_FULL:
            times, depth, sv = load_sv_smoothed(row['sv_zarr_path'], float(freq), *SMOOTHING[smooth])
            if threshold > -110:
                sv = np.where(sv >= threshold, sv, np.nan)
            quad = hv.QuadMesh((times, depth, sv), kdims=['ping_time', 'depth_m'], vdims=['Sv'])
            title = f"{freq:g} kHz - cleaned Sv, full resolution" + \
                (f", smoothed {smooth.split(' (')[0].lower()}" if smooth != 'Off' else '')
        else:
            mv = mv_all[np.isclose(mv_all['frequency_khz'], freq)]
            if mv.empty:
                continue
            grid = mv.pivot(index='depth_m', columns='ping_time', values='mvbs_db')
            values = grid.values if threshold <= -110 else np.where(grid.values >= threshold, grid.values, np.nan)
            quad = hv.QuadMesh((grid.columns.values, grid.index.values, values),
                               kdims=['ping_time', 'depth_m'], vdims=['Sv'])
            title = f"{freq:g} kHz - MVBS grid"
        slider = freq_sliders.get(freq)
        clim = tuple(slider.value) if slider is not None else (-85, -40)
        # Echogram convention: surface (0 m) at the top, depth increasing downward.
        img = rasterize(quad, aggregator='mean', dynamic=True).opts(
            cmap='viridis', clim=clim, colorbar=True, invert_yaxis=True, ylim=(0, ymax),
            responsive=True, height=height, active_tools=['pan'], title=title, tools=['hover', 'wheel_zoom'],
            xlabel='Ping time (UTC)' if k == len(freqs) - 1 else '', ylabel='Depth (m)',
            clabel='Sv (dB re 1 m-1)')
        plots.append(_decorate(img, row['file_key'], freq, ymax, **deco))
    if not plots:
        return pn.pane.Alert('No data for this selection.', alert_type='warning')
    return _stack(plots)


def _difference_impl(row, ymax, height, diff_range, deco):
    """Sv(high frequency) - Sv(low frequency) on the MVBS grid, where both have data."""
    freqs = _file_freqs(row)
    if len(freqs) < 2:
        return pn.pane.Alert('This recording has only one frequency.', alert_type='warning')
    f_hi, f_lo = max(freqs), min(freqs)
    mv = query_df("SELECT frequency_khz, ping_time, depth_m, mvbs_db FROM echogram_mvbs WHERE file_key = ?", [row['file_key']])
    diff = analysis.frequency_difference(mv, f_hi, f_lo)
    if diff.empty:
        return pn.pane.Alert('No cells have data at both frequencies after cleaning.', alert_type='warning')
    grid = diff.pivot(index='depth_m', columns='ping_time', values='delta_db')
    quad = hv.QuadMesh((grid.columns.values, grid.index.values, grid.values), kdims=['ping_time', 'depth_m'], vdims=['delta'])
    img = rasterize(quad, aggregator='mean', dynamic=True).opts(
        cmap='RdBu_r', clim=tuple(diff_range), colorbar=True, invert_yaxis=True, ylim=(0, ymax), responsive=True,
        height=height, active_tools=['pan'], tools=['hover', 'wheel_zoom'], xlabel='Ping time (UTC)', ylabel='Depth (m)',
        title=f"{f_hi:g} kHz minus {f_lo:g} kHz, MVBS grid (red = stronger at {f_hi:g} kHz)",
        clabel=f"Sv({f_hi:g}) - Sv({f_lo:g}) (dB)")
    return _stack([_decorate(img, row['file_key'], f_lo, ymax, **deco)])


def _echogram_impl(file_key, view, ymax, show_seabed, height, smooth, threshold, diff_range, show_layers, layer_m,
                   show_targets, show_schools, _version):
    row = _file_row(file_key)
    if row is None:
        return pn.pane.Alert('No processed data yet. Run main.py for a cruise first.', alert_type='warning')
    deco = dict(show_seabed=show_seabed, show_layers=show_layers, layer_m=layer_m, show_targets=show_targets,
                show_schools=show_schools)
    if view == VIEW_RAW:
        return _raw_level_impl(row, ymax, height, deco)
    if view == VIEW_DIFF:
        return _difference_impl(row, ymax, height, diff_range, deco)
    return _calibrated_impl(row, view, ymax, height, smooth, threshold, deco)


def _echogram_note(view, file_key, smooth):
    row = _file_row(file_key)
    found = ''
    if row is not None:
        found = f"  \nRecording **{row['filename']}** - frequencies detected: " + \
                ", ".join(f"{f:g} kHz" for f in _file_freqs(row)) + "."
    if view == VIEW_RAW:
        text = ("**Uncalibrated** echo level exactly as recorded: nothing removed, so the transmit pulse, "
                "noise and the area below the seabed all show. Each panel has its own colour scale; "
                "values are comparable within a panel only.")
    elif view == VIEW_DIFF:
        text = ("**Frequency difference** on the MVBS grid: the higher frequency minus the lower one, only where "
                "both have data above the noise. Red = stronger at the higher frequency. Interpret with care: the "
                "200 kHz channel is mostly noise below about 10 m in the Baja recording, so only shallow cells appear.")
    else:
        text = "Colour scale is **calibrated Sv**. Blank = removed (noise, near field, below seabed)."
        if view == VIEW_FULL and smooth != 'Off':
            text += (" **Smoothed:** values are averaged in blocks (power, not dB). A block with less than half "
                     "its samples left after cleaning stays blank, so where most data was judged to be noise "
                     "(for example 200 kHz below about 10 m here) the smoothed view is mostly blank.")
    return text + found


view_echogram = pn.Column(
    pn.bind(_echogram_note, view_select, file_select, smooth_select),
    pn.Card(details_pane, title='Cursor details (click an echogram)', collapsed=False, sizing_mode='stretch_width'),
    pn.bind(_echogram_impl, file_select, view_select, depth_max, seabed_toggle, plot_height, smooth_select,
            sv_threshold, diff_clim, layers_toggle, layer_select, targets_toggle, schools_toggle, clim_version),
    sizing_mode='stretch_width')


# 5. SEABED VIEW
def _seabed_impl(cruise):
    sb = query_df("SELECT file_key, ping_time, seabed_depth_m FROM echogram_seabed "
                  "WHERE cruise_id = ? ORDER BY ping_time", [cruise])
    if sb.empty:
        return pn.pane.Alert('No seabed detections for this cruise.', alert_type='warning')
    curve = hv.Curve(sb, kdims=['ping_time'], vdims=['seabed_depth_m']).opts(
        invert_yaxis=True, responsive=True, min_height=420, color='#00f2ff', line_width=1.5,
        xlabel='Ping time (UTC)', ylabel='Seabed depth (m)', title='Detected seabed depth',
        tools=['hover'], hooks=[_white_background])
    return pn.pane.HoloViews(curve, sizing_mode='stretch_both')


view_seabed = pn.Column(pn.bind(_seabed_impl, cruise_select), sizing_mode='stretch_both')


# 6. MAP OF SURVEY POSITIONS
def _to_mercator(lon, lat):
    """WGS84 degrees -> Web Mercator metres (the projection of the web map tiles)."""
    r = 6378137.0
    return np.radians(lon) * r, np.log(np.tan(np.pi / 4 + np.radians(lat) / 2)) * r


STATUS_COLORS = {'selected recording': '#ff1a1a', 'processed': '#ffb000', 'not yet processed': '#9aa0a6'}


def _limit_map_zoom(plot, element):
    """Ocean Base has tiles only to about zoom 10 in many regions; stop there instead of showing blank tiles."""
    for r in plot.state.renderers:
        src = getattr(r, 'tile_source', None)
        if src is not None:
            src.max_zoom = 10


def _map_impl(cruise, file_key):
    st = query_df("SELECT survey_id, latitude, longitude, position_confidence, processed "
                  "FROM echogram_stations WHERE cruise_id = ? ORDER BY survey_id", [cruise])
    st = st.dropna(subset=['latitude', 'longitude'])
    if st.empty:
        return pn.pane.Alert('No survey positions. Add cruises/<cruise>/stations.csv and rebuild.',
                             alert_type='warning')
    row = _file_row(file_key)
    selected = row['survey_id'] if row is not None and 'survey_id' in row else None
    st['x'], st['y'] = _to_mercator(st['longitude'].values, st['latitude'].values)
    st['status'] = np.where(st['survey_id'] == selected, 'selected recording',
                            np.where(st['processed'], 'processed', 'not yet processed'))
    st['position'] = st['position_confidence'].replace('', 'unrated')
    # Draw the selected recording last so it is never hidden behind a neighbouring dot.
    order = {'not yet processed': 0, 'processed': 1, 'selected recording': 2}
    st = st.sort_values('status', key=lambda c: c.map(order)).reset_index(drop=True)
    pts = hv.Points(st, kdims=['x', 'y'], vdims=['survey_id', 'latitude', 'longitude', 'position', 'status']).opts(
        color='status', cmap=STATUS_COLORS, size=13, line_color='white', line_width=1.5,
        tools=['hover'], show_legend=True, legend_position='bottom_left')
    # Esri Ocean Base (seafloor relief) with its labels; holoviews has no ready-made source for these two.
    ocean = (Tiles('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{Z}/{Y}/{X}',
                   name='EsriOceanBase') *
             Tiles('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Reference/MapServer/tile/{Z}/{Y}/{X}',
                   name='EsriOceanReference'))
    base = ocean * pts
    trk = query_safe("SELECT lat, lon FROM echogram_track WHERE cruise_id = ? AND lat IS NOT NULL ORDER BY ping_time", [cruise])
    if not trk.empty:
        tx, ty = _to_mercator(trk['lon'].values[::20], trk['lat'].values[::20])
        base = base * hv.Curve((tx, ty)).opts(color='#00f2ff', line_width=2)
    # Start with a sensible view.  A single position would otherwise zoom in so far that the
    # map has no tiles ("Map data not yet available"); show at least about 80 km around the points,
    # which is about the finest zoom (10) at which Ocean Base has tiles everywhere.
    # Centre on the stations; the track is still drawn but does not pull the view toward offshore legs.
    xs, ys = st['x'].to_numpy(), st['y'].to_numpy()
    half_min = 80000.0 / np.cos(np.radians(float(st['latitude'].mean())))      # Web Mercator metres for 80 km of ground
    cx, cy = (xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2
    hx, hy = max((xs.max() - xs.min()) * 0.6, half_min), max((ys.max() - ys.min()) * 0.6, half_min)
    plot = base.opts(
        responsive=True, min_height=620, xaxis=None, yaxis=None, active_tools=['wheel_zoom'],
        xlim=(cx - hx, cx + hx), ylim=(cy - hy, cy + hy),
        bgcolor='#0b1a2b',   # dark backdrop, shown if the map tiles cannot load
        title='Survey positions (source and accuracy of each position: see stations.csv)',
        hooks=[_white_border, _limit_map_zoom])
    return pn.pane.HoloViews(plot, sizing_mode='stretch_both')


view_map = pn.Column(pn.bind(_map_impl, cruise_select, file_select), sizing_mode='stretch_both')


# 7b. CSV EXPORT HELPER
def _csv_button(label, get_df, filename):
    """A download button that writes the table currently shown (``get_df()``) as CSV."""
    def _callback():
        df = get_df()
        buf = io.StringIO()
        (df if df is not None else pd.DataFrame()).to_csv(buf, index=False)
        buf.seek(0)
        return io.BytesIO(buf.getvalue().encode('utf-8'))
    return pn.widgets.FileDownload(callback=_callback, filename=filename, label=label, button_type='primary',
                                   embed=False, auto=True)


def _table(df, **kw):
    return pn.widgets.Tabulator(df, theme='midnight', disabled=True, sizing_mode='stretch_width', show_index=False,
                                pagination='remote', page_size=kw.pop('page_size', 12), **kw)


# 7. TABLES
def _mvbs_table_impl(file_key):
    df = query_df("SELECT ping_time, frequency_khz, depth_m, mvbs_db FROM echogram_mvbs "
                  "WHERE file_key = ? ORDER BY frequency_khz, ping_time, depth_m", [file_key])
    if df.empty:
        return pn.pane.Alert('No MVBS cells for the current selection.')
    return pn.Column(
        _csv_button('Download MVBS cells (CSV)', lambda: df, f"mvbs_{file_key}.csv"),
        pn.widgets.Tabulator(df, pagination='remote', page_size=15, theme='midnight',
                             header_filters=True, sizing_mode='stretch_both', disabled=True),
        sizing_mode='stretch_both')


view_mvbs_table = pn.Column(pn.bind(_mvbs_table_impl, file_select), sizing_mode='stretch_both')


def _summary_impl(cruise):
    files = query_df("SELECT filename, survey_id, latitude, longitude, first_ping, last_ping, n_pings, "
                     "frequencies_khz, waveform_mode, sound_speed_ms, seabed_min_m, seabed_max_m "
                     "FROM echogram_files "
                     "WHERE cruise_id = ? ORDER BY first_ping", [cruise])
    meta = query_df("SELECT * FROM build_metadata WHERE cruise_id = ?", [cruise])
    if files.empty:
        return pn.pane.Alert('No processed files for this cruise.')
    return pn.Column(
        pn.pane.Markdown("### Recordings"),
        pn.widgets.Tabulator(files, theme='midnight', disabled=True, sizing_mode='stretch_width'),
        pn.pane.Markdown("### Build provenance"),
        pn.widgets.Tabulator(meta.T.reset_index().rename(columns={'index': 'field', 0: 'value'}),
                             theme='midnight', disabled=True, sizing_mode='stretch_width',
                             show_index=False),
        sizing_mode='stretch_both')


view_summary = pn.Column(pn.bind(_summary_impl, cruise_select), sizing_mode='stretch_both')



# 7c. REGIONS (draw boxes, get statistics, export as annotations)
region_classes = {}                        # region number -> user-entered class name


def _regions_impl(file_key):
    row = _file_row(file_key)
    if row is None:
        return pn.pane.Alert('No processed data yet.', alert_type='warning')
    mv = query_df("SELECT frequency_khz, ping_time, depth_m, mvbs_db FROM echogram_mvbs WHERE file_key = ?", [file_key])
    if mv.empty:
        return pn.pane.Alert('No MVBS cells for this recording.', alert_type='warning')
    f_lo = float(mv['frequency_khz'].min())
    g = mv[np.isclose(mv['frequency_khz'], f_lo)].pivot(index='depth_m', columns='ping_time', values='mvbs_db')
    vmin, vmax = np.nanpercentile(g.values, [1, 99.5])
    quad = hv.QuadMesh((g.columns.values, g.index.values, g.values), kdims=['ping_time', 'depth_m'], vdims=['Sv'])
    img = rasterize(quad, aggregator='mean', dynamic=True).opts(
        cmap='viridis', clim=(float(vmin), float(vmax)), colorbar=True, invert_yaxis=True, ylim=(0, depth_max.value),
        responsive=True, height=420, tools=['wheel_zoom', 'pan'], title=f"{f_lo:g} kHz MVBS - drag with the box tool to draw regions",
        xlabel='Ping time (UTC)', ylabel='Depth (m)')
    t0 = g.columns.values[0]
    rects = hv.Rectangles([(t0, 0.0, t0, 0.0)], kdims=['x0', 'y0', 'x1', 'y1']).opts(
        fill_alpha=0.15, fill_color='white', line_color='white', line_width=2)
    edit = hv.streams.BoxEdit(source=rects, num_objects=10, styles={})
    table = pn.widgets.Tabulator(pd.DataFrame(), theme='midnight', sizing_mode='stretch_width', show_index=False,
                                 editors={'frequency_khz': None, 'n_cells': None, 'mean_sv_db': None, 'sa_m2_m2': None,
                                          'nasc_m2_nmi2': None, 'delta_hi_minus_lo_db': None, 'start_utc': None,
                                          'end_utc': None, 'top_m': None, 'bottom_m': None, 'region': None, 'class': 'input'})
    msg = pn.pane.Markdown('Hold **Shift** and drag on the echogram to draw a region (the box-edit tool is already active; click a box and press Backspace to delete it, drag to move it). '
                           'Type a class name (for example "krill", "fish", "noise") in the table. '
                           'Statistics use the MVBS grid at every frequency.')

    def refresh(data):
        data = data or {}
        rows = []
        n = 0                                   # region number (the empty placeholder box is skipped)
        for (x0, x1, y0, y1) in (zip(data.get('x0', []), data.get('x1', []), data.get('y0', []), data.get('y1', []))):
            ta, tb = _to_datetime64(x0), _to_datetime64(x1)
            if ta == tb or y0 == y1:
                continue
            n += 1
            st = analysis.region_stats(mv, ta, tb, y0, y1)
            for r in st.to_dict('records'):
                rows.append(dict(region=n, start_utc=str(pd.Timestamp(min(ta, tb)).round('s')),
                                 end_utc=str(pd.Timestamp(max(ta, tb)).round('s')),
                                 top_m=round(min(y0, y1), 1), bottom_m=round(max(y0, y1), 1),
                                 **{k: (float(f'{v:.4g}') if isinstance(v, float) else v) for k, v in r.items()},
                                 **{'class': region_classes.get(n, '')}))
        table.value = pd.DataFrame(rows) if rows else pd.DataFrame(
            columns=['region', 'start_utc', 'end_utc', 'top_m', 'bottom_m', 'frequency_khz', 'n_cells', 'mean_sv_db',
                     'sa_m2_m2', 'nasc_m2_nmi2', 'class'])

    def classes_changed(event):
        for r in table.value.to_dict('records'):
            region_classes[r['region']] = r['class']

    edit.add_subscriber(refresh)
    table.param.watch(classes_changed, 'value')
    plot = (img * rects).opts(active_tools=['box_edit'], hooks=[_white_background])
    return pn.Column(msg, pn.pane.HoloViews(plot, sizing_mode='stretch_width'), table,
                     _csv_button('Download regions (CSV)', lambda: table.value, f"regions_{file_key}.csv"),
                     sizing_mode='stretch_width')


view_regions = pn.Column(pn.bind(_regions_impl, file_select), sizing_mode='stretch_width')


# 7d. INTEGRATION (layers, sa, NASC, echo counts)
int_layer = pn.widgets.Select(name='Layer thickness (m)', options=[5, 10, 20, 50], value=10, width=160)
int_time = pn.widgets.Select(name='Time bin', options={'Whole recording': 'all', 'MVBS time step (10 s)': '10s', '1 min': '1min', '5 min': '5min'},
                             value='all', width=160)


def _integration_df(file_key, layer_m, time_bin):
    mv = query_df("SELECT frequency_khz, ping_time, depth_m, mvbs_db FROM echogram_mvbs WHERE file_key = ?", [file_key])
    whole = time_bin == 'all'
    out = analysis.integrate_layers(mv, float(layer_m), '3650D' if whole else time_bin)
    if out.empty:
        return out
    out['time_bin'] = 'whole recording' if whole else out['time_bin'].astype(str)
    for c in ('coverage', 'mean_sv_db', 'nasc_m2_nmi2'):
        out[c] = out[c].astype(float).round(4)
    out['sa_m2_m2'] = out['sa_m2_m2'].astype(float).map(lambda v: f'{v:.4g}')   # text, so tiny values keep their digits
    tg = query_safe("SELECT frequency_khz, depth_m FROM echogram_targets WHERE file_key = ?", [file_key])
    if not tg.empty:
        tg['layer_top_m'] = np.floor(tg['depth_m'] / float(layer_m)) * float(layer_m)
        cnt = tg.groupby(['frequency_khz', 'layer_top_m']).size().rename('echo_count').reset_index()
        out = out.merge(cnt, how='left', on=['frequency_khz', 'layer_top_m'])
        # counts are per recording, so they are not repeated for each time bin
        out['echo_count'] = np.nan if not whole else out['echo_count'].fillna(0).astype(int)
    return out


def _integration_impl(file_key, layer_m, time_bin):
    df = _integration_df(file_key, layer_m, time_bin)
    if df.empty:
        return pn.pane.Alert('No MVBS cells for this recording.', alert_type='warning')
    prof = df.groupby(['frequency_khz', 'layer_top_m'], as_index=False).agg(nasc=('nasc_m2_nmi2', 'mean'))
    curves = []
    for f, g in prof.groupby('frequency_khz'):
        curves.append(hv.Curve(g, 'nasc', 'layer_top_m', label=f"{f:g} kHz").opts(line_width=2, tools=['hover']))
    plot = hv.Overlay(curves).opts(invert_yaxis=True, logx=True, responsive=True, height=340, xlabel='NASC (m2 / nmi2), log scale',
                                   ylabel='Layer top (m)', title=f"Mean NASC by {layer_m} m layer", legend_position='right',
                                   hooks=[_white_background])
    return pn.Column(pn.pane.HoloViews(plot, sizing_mode='stretch_width'), _table(df, page_size=15),
                     _csv_button('Download layers (CSV)', lambda: _integration_df(file_select.value, int_layer.value, int_time.value),
                                 f"layers_{file_key}.csv"), sizing_mode='stretch_width')


view_integration = pn.Column(
    pn.pane.Markdown("**Echo integration** of the cleaned MVBS grid into depth layers. sa is the integrated backscatter "
                     "(m2/m2); NASC = 4 pi 1852^2 sa (m2/nmi2). Coverage is the share of a layer's cells that held data after "
                     "cleaning - treat layers with low coverage with caution. Note: the 200 kHz channel in this recording is mostly noise below about 10 m, so its layer values are not trustworthy there. `echo_count` is the number of single targets "
                     "detected in the layer (whole-recording view only)."),
    pn.Row(int_layer, int_time),
    pn.bind(_integration_impl, file_select, int_layer, int_time), sizing_mode='stretch_width')


# 7e. SINGLE TARGETS AND SCHOOLS
def _targets_impl(file_key):
    t = query_safe("SELECT * EXCLUDE (cruise_id) FROM echogram_targets WHERE file_key = ? ORDER BY frequency_khz, ping_time", [file_key])
    sc = query_safe("SELECT * EXCLUDE (cruise_id) FROM echogram_schools WHERE file_key = ? ORDER BY frequency_khz, school_id", [file_key])
    parts = [pn.pane.Markdown("**Single targets** are isolated echoes that pass the thresholds in the config (TS, SNR, echo length, "
                              "beam compensation, peak above local background). These are *candidates*: without ground truth "
                              "(net or camera) they are not confirmed fish. TS is compensated for the beam pattern only where split-beam "
                              "angles exist (38 kHz here); `compensated` says which.")]
    if t.empty:
        parts.append(pn.pane.Alert('No single targets detected in this recording.', alert_type='info'))
    else:
        h = hv.Histogram(np.histogram(t['ts_comp_db'], bins=20)).opts(xlabel='TS (dB re 1 m2)', ylabel='Count', title='TS distribution',
                                                                      responsive=True, height=260, tools=['hover'], hooks=[_white_background])
        sc_plot = hv.Scatter(t, 'ts_comp_db', ['depth_m', 'frequency_khz']).opts(invert_yaxis=True, size=9, color='frequency_khz', cmap='Category10',
                                                                                  responsive=True, height=260, xlabel='TS (dB)', ylabel='Depth (m)',
                                                                                  title='TS versus depth', tools=['hover'], hooks=[_white_background])
        parts += [pn.Row(pn.pane.HoloViews(h, sizing_mode='stretch_width'), pn.pane.HoloViews(sc_plot, sizing_mode='stretch_width')),
                  _table(t), _csv_button('Download targets (CSV)', lambda: t, f"targets_{file_key}.csv")]
    parts.append(pn.pane.Markdown("**Detected schools / layers** - connected areas of smoothed backscatter above the threshold."))
    if sc.empty:
        parts.append(pn.pane.Alert('No schools or layers detected.', alert_type='info'))
    else:
        parts += [_table(sc), _csv_button('Download schools (CSV)', lambda: sc, f"schools_{file_key}.csv")]
    return pn.Column(*parts, sizing_mode='stretch_width')


view_targets = pn.Column(pn.bind(_targets_impl, file_select), sizing_mode='stretch_width')


def _unchecked_banner(cruise):
    """A prominent warning while the cruise's settings are untouched defaults (CONFIG_CHECKED = false)."""
    cfg = query_safe("SELECT config_json FROM echogram_config WHERE cruise_id = ?", [cruise])
    if cfg.empty:
        return None
    c = json.loads(cfg['config_json'].iloc[0])
    if not c.get('_CONFIG_UNCHECKED'):
        return None
    return pn.pane.Alert(
        "**Default settings - check them before you trust these results.** This cruise was built with default "
        "processing settings that were not tuned for its data. Check that the red **seabed line** follows the real "
        "seabed, that noisy channels and depths are removed (and real echoes are not), and that the **water values** "
        "below are measured, not the EK80 software defaults. Then edit the cruise's `config.toml`, rebuild, and set "
        "`CONFIG_CHECKED = true` to remove this warning.", alert_type='danger', sizing_mode='stretch_width')


# 7f. CALIBRATION AND QC
def _calibration_impl(file_key, _cruise=None):
    cal = query_safe("SELECT * EXCLUDE (cruise_id, file_key) FROM echogram_calibration WHERE file_key = ? ORDER BY frequency_khz", [file_key])
    qc = query_safe("SELECT * EXCLUDE (cruise_id, file_key) FROM echogram_qc WHERE file_key = ? ORDER BY frequency_khz", [file_key])
    if cal.empty and qc.empty:
        return pn.pane.Alert('No calibration or QC tables yet - rebuild the database with the current main.py.', alert_type='warning')
    parts = []
    banner = _unchecked_banner(cruise_select.value)
    if banner is not None:
        parts.append(banner)
    if not cal.empty:
        if (cal['environment_source'] != 'config.toml').any():
            parts.append(pn.pane.Alert('Sound speed, absorption and salinity come from built-in defaults, not from this cruise. '
                                       'Set [environment] in config.toml (for example from a CTD cast) for trustworthy absolute Sv.', alert_type='warning'))
        parts += [pn.pane.Markdown("### Calibration and environment applied"),
                  _table(cal.set_index('frequency_khz').T.reset_index().rename(columns={'index': 'parameter'}), page_size=30),
                  pn.pane.Markdown("Gain, Sa correction and equivalent beam angle are the values stored in the file by the EK80 software. "
                                   "No sphere calibration has been applied by this code. Transceiver impedance is only stored once per file, "
                                   "so a per-ping impedance trace cannot be produced here."),
                  _csv_button('Download calibration (CSV)', lambda: cal, f"calibration_{file_key}.csv")]
    if not qc.empty:
        parts += [pn.pane.Markdown("### Data quality"), _table(qc.set_index('frequency_khz').T.reset_index().rename(columns={'index': 'check'}), page_size=30),
                  _csv_button('Download QC (CSV)', lambda: qc, f"qc_{file_key}.csv")]
    return pn.Column(*parts, sizing_mode='stretch_width')


view_calibration = pn.Column(pn.bind(_calibration_impl, file_select, cruise_select), sizing_mode='stretch_width')


# 7g. PROCESSING (dataflow)
def _processing_impl(cruise):
    cfg = query_safe("SELECT config_json FROM echogram_config WHERE cruise_id = ?", [cruise])
    if cfg.empty:
        return pn.pane.Alert('No stored configuration for this cruise.', alert_type='warning')
    c = json.loads(cfg['config_json'].iloc[0])
    g = c.get
    steps = [
        ("1. Read", "Raw EK80 .raw files, converted with Echopype (an optional .nc archive can be written with --write-netcdf; it is never read back)."),
        ("2. Calibrate", f"Sv from power, gain and Sa correction; T {g('TEMPERATURE_C')} C, S {g('SALINITY_PSU')} PSU, "
                         f"{g('PRESSURE_DBAR')} dbar, sound speed {g('SOUND_SPEED_MS')} m/s."),
        ("3. Excluded pings", f"Time ranges blanked: {g('EXCLUDE_TIME_RANGES_UTC') or 'none'}."),
        ("4. Near field", f"Remove the top {g('NEAR_FIELD_EXCLUSION_M')} m."),
        ("5. Seabed", f"Detected on {str(g('SEABED_CHANNEL_KHZ')) + ' kHz' if g('SEABED_CHANNEL_KHZ') else 'each file’s lowest-frequency channel (automatic)'} (threshold {g('SEABED_THRESHOLD_DB')} dB, offset {g('SEABED_OFFSET_M')} m, "
                      f"margin {g('SEABED_MASK_MARGIN_M')} m); everything below is blank."),
        ("6. Background noise", f"{'On' if g('BACKGROUND_NOISE') else 'Off'}: {g('BN_PING_NUM')} pings x {g('BN_RANGE_M')} m, "
                                f"SNR threshold {g('SNR_THRESHOLD_DB')} dB."),
        ("7. Impulse noise", f"{'On' if g('IMPULSE_NOISE') else 'Off'}: {g('IMPULSE_THRESHOLD_DB')} dB over neighbours "
                             f"(+/-{g('IMPULSE_SIDE_PINGS')} pings)."),
        ("8. Transient noise", f"{'On' if g('TRANSIENT_NOISE') else 'Off (default)'}: {g('TRANSIENT_THRESHOLD_DB')} dB over local median "
                               f"(+/-{g('TRANSIENT_SIDE_PINGS')} pings, ignoring top {g('TRANSIENT_EXCLUDE_ABOVE_M')} m)."),
        ("9. Grid", f"MVBS cells {g('MVBS_PING_TIME_BIN')} x {g('MVBS_RANGE_BIN_M')} m (linear averaging)."),
        ("10. Detect", f"Single targets {'on' if g('TARGET_DETECTION') else 'off'}; schools {'on' if g('SCHOOL_DETECTION') else 'off'} "
                       f"({g('SCHOOL_THRESHOLD_DB')} dB, {g('SCHOOL_MIN_THICKNESS_M')} m x {g('SCHOOL_MIN_LENGTH_S')} s)."),
    ]
    src = g('_CONFIG_SOURCE') or 'not recorded (built before this was tracked)'
    banner = _unchecked_banner(cruise) or pn.pane.Markdown(f"Configuration used: **{src}**.")
    flow = pn.pane.Markdown('\n\n'.join(f"**{a}** - {b}" for a, b in steps))
    flat = pd.DataFrame(sorted((k, str(v)) for k, v in c.items()), columns=['setting', 'value'])
    return pn.Column(pn.pane.Markdown("### Processing steps, in the order they were applied"), banner, flow,
                     pn.pane.Markdown("### All settings stored with this build"), _table(flat, page_size=20), sizing_mode='stretch_width')


view_processing = pn.Column(pn.bind(_processing_impl, cruise_select), sizing_mode='stretch_width')


# 8. ASSEMBLY
tabs = pn.Tabs(
    ("Echogram", view_echogram),
    ("Regions", view_regions),
    ("Integration", view_integration),
    ("Targets & Schools", view_targets),
    ("Seabed Depth", view_seabed),
    ("Map", view_map),
    ("Calibration & QC", view_calibration),
    ("Processing", view_processing),
    ("MVBS Table", view_mvbs_table),
    ("Cruise Summary", view_summary),
    dynamic=True, active=0)

dashboard = pn.template.FastListTemplate(
    title="R/V Western Flyer - Echogram (EK80)",
    sidebar=[cruise_select, file_select, pn.pane.Markdown("**View**"), view_select, smooth_select, sliders_box,
             lock_clim, sv_threshold, diff_clim, depth_max, plot_height,
             reset_button, pn.pane.Markdown("**Overlays**"), seabed_toggle, schools_toggle, targets_toggle, layers_toggle, layer_select],
    main=[tabs], accent_base_color="#00f2ff", header_background="#1a1a1a")

dashboard.servable()
