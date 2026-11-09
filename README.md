# Western Flyer EK80 Echograms

This repository contains the processing pipeline and visualization suite for Simrad EK80 echosounder data collected aboard the R/V Western Flyer.

## Feature Summary

* **Raw `.raw` Processing:** The pipeline starts from the instrument's native `.raw` files and uses [Echopype](https://echopype.readthedocs.io) for conversion and calibration. No commercial echosounder software is needed.
* **Optional lossless archive:** With `--write-netcdf`, each `.raw` file is also repackaged, with no processing, as a SONAR-netCDF4 file (`.nc`, ICES convention). It is off by default because the files are large and the dashboard never reads them.
* **Calibrated Sv:** Complex backscatter is converted to volume backscattering strength (Sv, dB re 1 m-1) using the gain, equivalent beam angle and sa-correction stored in the file, spreading loss and sound absorption. The calibration mode (CW or broadband) is chosen from what was actually transmitted.
* **Standard cleaning chain (per channel):** near-field exclusion, then seabed detection, then background-noise removal (De Robertis & Higginbottom 2007), then impulse-noise masking (Ryan et al. 2015), then masking below the seabed.
* **Two data products:** cleaned full-resolution Sv (Zarr, for zooming) and an MVBS grid (mean volume backscattering strength per range-bin x ping-time-bin cell, for analysis).
* **Modular Cruise Architecture:** Each expedition is isolated under `cruises/<cruise_id>/`. Adding a cruise requires only a directory with the `.raw` files and a `config.toml`; the directory name is the cruise ID.
* **Atomic, Idempotent Storage:** Small tables go to DuckDB (and, with `--write-parquet`, to Parquet). Per-cruise records are replaced atomically inside a single transaction. A `build_metadata` table records the SHA-256 of `config.toml` and of every raw file, the Echopype version and the git commit, so every row is traceable.
* **Per-Build Audit Logging:** Each invocation writes a fresh, timestamped log file (`logs/wf_echogram_build_<cruise>_<timestamp>.log`).
* **Survey positions:** The `.raw` files carry no GPS, so positions come from `stations.csv` (one row per survey) and are joined onto each file by `station_id`.
* **Echosounder analysis tools (all in Python):** frequency difference, echo integration by depth layer (sa, NASC), single-target detection (TS, split-beam angles), school/layer detection, drawn regions with statistics and a class label, excluded-ping and transient-noise handling, calibration and QC tables, an optional ship track from a navigation log, and CSV export of every table.
* **Interactive Dashboard:** A Panel/HoloViews dashboard renders full-resolution echograms with Datashader (fast at any zoom), the MVBS grid, the seabed line, a map of survey positions, tables and build provenance.

## Root Directory Notes

The project root can be named anything (e.g., `westernflyer_cruise_echogram`).

To avoid Windows permission issues, place the root folder inside your Windows user directory:
`C:\Users\<your_username>\westernflyer_cruise_echogram`

## Python Installation Requirements

This project requires **Python 3.13 or newer** (Echopype 0.12 does not install on 3.12).

### Windows

Install Python 3.13 from python.org (not from the Microsoft Store, which puts it in a different folder than the instructions below assume).

1. If you installed Python 3.13 from the Microsoft Store, remove it first: Windows Settings, Apps, "Python 3.13", Uninstall.
2. Download the "Windows installer (64-bit)" for Python 3.13 from https://www.python.org/downloads/windows/.
3. Run it and choose **Install Now**. Leave **"Install for all users"** unticked, so it installs for your user only. Leave the **py launcher** ticked.
4. When the installer asks to **"Add python.exe to PATH"**, select **No** (leave it unticked).
5. Python installs to `C:\Users\<your_username>\AppData\Local\Programs\Python\Python313\`.
6. Check it in Command Prompt: `py -0p` should list `Python313` under `AppData\Local\Programs\Python`.

### macOS (Homebrew)

    brew install python@3.13

### Linux (Ubuntu / Debian) using apt install

    sudo add-apt-repository ppa:deadsnakes/ppa
    sudo apt update
    sudo apt install python3.13 python3.13-venv

### Linux (Ubuntu / Debian) using pyenv

This has the advantage of not requiring root privileges, and it does not install Python 3.13 system wide.

    # Install pyenv:
    curl -fsSL https://pyenv.run | bash
    # Install Python 3.13:
    pyenv install 3.13

> The pipeline was developed on Linux with Python 3.13 and Echopype 0.12.0. The 20250416_Baja build has also been run successfully on Windows with Python 3.13 and Echopype 0.12.0. The dashboard has been tested on Linux; it has not yet been confirmed on Windows.

## Project Structure

*`logs/` and `processed/` are created automatically on first run.*

    .
    |   echogram_holoviews.py   dashboard
    |   echogram_processing.py  noise removal, seabed, masking, MVBS gridding, QC, schools
    |   echogram_analysis.py    layers, NASC, regions, difference, school and track helpers
    |   echogram_targets.py     single-target (TS) detection
    |   ek80_ingestion.py       .raw -> EchoData (optional .nc), calibration to Sv
    |   main.py                 orchestrator
    |   pyproject.toml
    +---cruises
    |   \---20250416_Baja
    |       |   config.toml     processing parameters
    |       |   stations.csv    survey positions (station_id, lat, lon)
    |       \---raw
    |               250416WF-D20250416-T191232 (1).raw
    +---logs
    +---processed
    |   |   wf_echogram.duckdb
    |   \---20250416_Baja
    |       |   mvbs.parquet, seabed.parquet, files.parquet, stations.parquet   (only with --write-parquet)
    |       +---netcdf          <file_key>.nc    lossless archive of the raw file (only with --write-netcdf)
    |       +---sv              <file_key>.zarr  cleaned full-resolution Sv
    |       \---raw_level       <file_key>.zarr  uncalibrated echo level, all frequencies
    \---scripts
        +---linux_mac   echogram_build.sh, echogram_build_all.sh, echogram_build_20250416_Baja.sh, echogram_build_Tutorial_IntroToLiveViewing.sh, echogram_dashboard.sh
        \---windows     echogram_build.bat, echogram_build_all.bat, echogram_build_20250416_Baja.bat, echogram_build_Tutorial_IntroToLiveViewing.bat, echogram_dashboard.bat

## Cruises

| Cruise ID | Dates | Raw files | Notes |
|---|---|---|---|
| `20250416_Baja` | Apr 16 2025 (recording starts 19:12 UTC) | 1 | 38 and 200 kHz, CW pings, 389 pings, seabed 84-89 m. The `.raw` file has no valid GPS position; positions come from `stations.csv`. |
| `Tutorial_IntroToLiveViewing` | May 13 2015 (file names start 09:09:35 to 09:11:05, about 2 minutes) | 4 | Echoview "Getting Started" tutorial data. 120 kHz (ES120-7C, 250 W) and 200 kHz (ES200-7C, 150 W), split-beam, **sequential pinging** (the channels transmit in turn) with **broadband FM chirps** of 2.048 ms (95-160 and 160-260 kHz). Shallow: seabed about 18 m in the part tested. The EK80 water values in the files are the factory defaults (10 C, 35 PSU, 100 m). The files do contain GPS: the ship was stationary at about 59.446 N, 10.380 E (Oslofjord, Norway), and `stations.csv` was made from those fixes. Its `config.toml` holds default settings and is not yet checked (`CONFIG_CHECKED = false`), so the dashboard shows a warning banner until it is. |

## File Descriptions

* **main.py** - Orchestrator. Parses CLI args, runs the per-file pipeline, writes the Zarr and DuckDB outputs (plus the optional `.nc` archive and Parquet copies), and stamps build provenance.
* **ek80_ingestion.py** - Ingestion and calibration. Reads `.raw` via Echopype, can write the optional `.nc` archive, applies environment overrides and calibrates to Sv. Also loads `stations.csv`.
* **echogram_processing.py** - Cleaning and gridding. Noise removal, seabed detection and masking, MVBS gridding, flattening results to tables.
* **echogram_analysis.py** - Pure NumPy/pandas helpers: block averaging, layer integration (sa, NASC), region statistics, frequency difference, school detection, navigation interpolation, transient-noise mask.
* **echogram_targets.py** - Single-target detection from TS, Sv and split-beam angles.
* **echogram_holoviews.py** - The Panel dashboard.
* **scripts/** - Windows `.bat` and Linux/macOS `.sh` wrappers for building and for launching the dashboard.

## Setup: Virtual Environment and Installation

### Windows (Command Prompt)

    cd C:\Users\<your_username>\westernflyer_cruise_echogram
    C:\Users\<your_username>\AppData\Local\Programs\Python\Python313\python.exe -m venv .venv
    .venv\Scripts\activate
    pip install -e .

### macOS / Linux

    cd westernflyer_cruise_echogram
    python3.13 -m venv .venv
    source .venv/bin/activate
    pip install -e .

## Configuration: defaults and automatic config.toml

Settings live in `cruises/<cruise_id>/config.toml`. You do not have to write one:

* **A cruise with no `config.toml`** gets one written automatically at the start of the build, with every setting at its default, from the first `.raw` file: the channels and pulse types found, the water values stored in the file (flagged when they look like the EK80 factory defaults of 10 C, 35 PSU, 100 m), and the seabed and school channel set to the lowest frequency. The file contains `CONFIG_CHECKED = false`. You can also write it without building: `python main.py <cruise_id> --init-config`. An existing file is never overwritten.
* **Defaults must be checked.** While `CONFIG_CHECKED = false`, the build prints a **REVIEW REQUIRED** warning at the start and end, and the dashboard's *Calibration & QC* and *Processing* tabs show a red warning banner. The defaults were tuned on one recording, so check that the red seabed line follows the real seabed, that noisy channels and depths are removed without removing real echoes, and that the water values are measured. Edit the file, rebuild, and set `CONFIG_CHECKED = true` to stop the warnings. A config without that key (such as the Baja one) is treated as checked.
* **`config.defaults.toml`** (project folder) holds the default for every setting. A key a cruise's `config.toml` leaves out falls back to it, so a cruise file can be just the few lines you change.

## Preparing a New Cruise

1. Create `cruises/<cruise_id>/raw/` and put the `.raw` files in it.
2. Skip this step to have the build write a default-settings `config.toml` for you (then check it; see *Configuration* above), run `python main.py <cruise_id> --init-config` to write it now, or copy one from an existing cruise.
3. Optionally add `nav.csv` (columns `time_utc`, `lat`, `lon`, decimal degrees) to get a per-ping ship track, distance and speed. Add `stations.csv` (see below) with one row per survey. Without it the build still works, but positions are empty and a warning is logged.
4. **Set the `[environment]` values** (temperature, salinity, pressure) from your CTD data. See the notes under *Known limitations*.
5. Run the build (next section).

## stations.csv

One row per survey. `station_id` is the start of the raw file name, up to the first `-` (for `250416WF-D20250416-T191232 (1).raw` it is `250416WF`).

| Column | Meaning |
|---|---|
| `station_id` | Required, unique |
| `lat`, `lon` | Required. Signed decimal degrees (north positive, west negative; blank if unknown). The same columns and format as the CTD project's `stations.csv` |
| `confidence`, `note` | Optional, carried into the database and shown on the map |

**stations.csv is written automatically from GPS when it is missing.** If a cruise has no `stations.csv`, the build reads the GPS fixes recorded inside the raw files (NMEA `GGA` sentences in the navigation datagrams) and writes one row per station id: the median fix, a confidence from the fix quality (`good`, or `fair` for plain GPS or a moving ship), and a note with the fix count and spread. It never overwrites a file you made, and it logs a warning asking you to check the positions. If the ship moved more than 500 m during a survey, one position is a poor description, so the build also writes a time-based `nav.csv` (one fix per second, timed with the same clock as the pings, so a GPS clock offset does not matter) and the dashboard draws the track. If the files hold no usable GPS (the Baja file does not), positions stay empty and the build says so.

**Positions from GPS in the raw files.** Some `.raw` files carry the ship's GPS as NMEA text (`$..GGA` sentences). `Tutorial_IntroToLiveViewing` does: its `stations.csv` row is the median of the 243 fixes found in the four files (all fix quality 4, the ship did not move: about 9 m spread in latitude, under 1 m in longitude), rated `good`. All four files share the station id `EK80_SimradEcho_WC381_Sequential`, so one row covers them. The GPS clock in those files reads about 5 minutes later than the file names, which matters only if you match ping times to GPS times. A moving survey needs a time-based `nav.csv` (`time_utc`, `lat`, `lon`) instead of one position per survey; the build writes one automatically (see above).

The 20250416_Baja positions were read off the "Golfo de California 2025 - Echosounder surveys" map (a picture with no coordinates), so they are approximate (about +/- 3 km). Replace them with positions from the ship's navigation log when you have them.

## Running the Pipeline

### Windows (Command Prompt)

    scripts\windows\echogram_build.bat 20250416_Baja
    scripts\windows\echogram_build_Tutorial_IntroToLiveViewing.bat
    scripts\windows\echogram_build_all.bat

### macOS / Linux

Make the scripts executable once (first time only):

    chmod +x scripts/linux_mac/*.sh

Then run:

    ./scripts/linux_mac/echogram_build.sh 20250416_Baja
    ./scripts/linux_mac/echogram_build_Tutorial_IntroToLiveViewing.sh

Or directly: `python main.py 20250416_Baja`. The Baja build takes about half a minute for the one 105 MB file.

### Messages you can ignore

During a build, Echopype prints these lines; they are notices, not errors:

* `EchoData object (converted raw file) does not contain valid Platform location data` and `xarray Dataset does not contain valid location data`: the `.raw` file has no GPS fixes, so Echopype cannot attach a location to its own output. Survey positions come from `stations.csv` instead and are added separately.
* `parsing file ...` and `saving ...`: progress messages.

A real problem appears as a line starting `WARNING` or `ERROR` in the build output and in `logs/`.

### CLI arguments

| Argument | Meaning |
|---|---|
| `cruise_id` | Folder name under `cruises/` |
| `--range-bin 1.0` | Override `MVBS_RANGE_BIN_M` (m) |
| `--ping-bin 10s` | Override `MVBS_PING_TIME_BIN` |
| `--raw-dir PATH` | Raw files elsewhere (default `cruises/<id>/raw/`) |
| `--db PATH` | DuckDB output path (default `processed/wf_echogram.duckdb`) |
| `--no-sv-store` | Skip the full-resolution Zarr stores (MVBS only; much smaller output) |
| `--write-netcdf` | Also save each `.raw` file as a lossless SONAR-netCDF4 archive (`processed/<cruise_id>/netcdf/<file_key>.nc`). Off by default; adds build time and about as much disk as the raw file |
| `--write-parquet` | Also write the tables as Parquet files (`processed/<cruise_id>/*.parquet`). Off by default |
| `--verbose` | DEBUG logging |

### Failed files

A file that fails to process is logged and skipped; the rest of the cruise still builds. The exit code is 1 only if every file failed.

## Running the Dashboard

    scripts\windows\echogram_dashboard.bat

or `python -m panel serve echogram_holoviews.py --port 5008 --show`. The dashboard runs on port 5008 (http://localhost:5008/echogram_holoviews), so it can run alongside the CTD dashboard (5006) and the flow-through dashboard (5007). Choose the cruise and the recording in the sidebar. The dashboard detects which frequencies the recording contains (38 and 200 kHz for the Baja file) and stacks one panel per frequency on a shared time axis. The *Map* tab shows every survey in `stations.csv` on an Esri Ocean Base map showing seafloor relief (needs internet), with the selected recording highlighted in red.

The Echogram tab has four views (the fourth, *Frequency difference (MVBS)*, is described below):

* **Full resolution (cleaned Sv):** every ping and range sample after cleaning, rasterized by Datashader. Calibrated Sv; blank where data was removed (near field, noise, spikes, below the seabed).
* **MVBS grid (gridded mean):** the same data averaged into cells of `MVBS_RANGE_BIN_M` by `MVBS_PING_TIME_BIN` (1 m by 10 s by default), so it is blocky by design. Good for tables and comparing recordings, not for looking at detail.
* **Uncalibrated echo level:** the data as recorded, every detected frequency stacked, with no gain, noise removal or seabed masking. Each panel has its own colour scale, so compare shapes within a panel, not levels between them. Use it for a first look; use the cleaned views for anything quantitative.

**Smoothing (full resolution view only).** The cleaned full-resolution view draws every sample, so it looks grainy and has white holes where noise was removed. The **Smoothing** selector averages the cleaned data in blocks (Light 0.25 m x 3 pings, Medium 0.5 m x 5 pings, Strong 1 m x 9 pings) before drawing, which makes layers easier to follow. Averaging is done on power rather than dB, and a block with fewer than half of its samples left after cleaning stays blank, so smoothing never fills in areas that were mostly removed. Where most of a channel was judged to be noise (200 kHz below about 10 m in the Baja recording, where only about 10 % of samples survive), the smoothed view is mostly blank. Smoothing only changes the display; the stored data is unchanged. The colour ranges are refitted when you change it.

**Making the plots bigger.** Use **Plot height per frequency (px)** in the sidebar (250 to 1000); the page scrolls when the panels are taller than the window. The sidebar can be hidden with the menu button at the top left, which gives the plots the full width. The mouse wheel scrolls the page over a plot (the default plot tool is pan); click the wheel-zoom icon in the plot toolbar if you want the wheel to zoom instead. **Maximum depth shown** sets the depth range.

**Colour ranges.** There is one colour-range slider per detected frequency, and each is set automatically from the data (1st to 99.5th percentile) whenever you change the recording or the view, because 38 and 200 kHz have very different Sv levels. You can still move them by hand. They do not apply to the uncalibrated view.

Tick **Lock colour ranges (compare recordings)** to stop the automatic reset. The dashboard then remembers one range per frequency (set the sliders the way you want them and each comes back whenever you select a recording) and uses them for every recording and for both calibrated views. That makes recordings directly comparable. Unticking the box returns to automatic ranges.

For long recordings, use the MVBS view or zoom into a time window. The full-resolution view loads every frequency of the selected file into memory (about 33 MB per frequency for 83 s at 200 kHz).

**Cursor details.** Click any echogram: a cyan cross-hair appears on every frequency and the *Cursor details* card shows the time, depth, ping number, cleaned Sv, uncalibrated level, seabed depth, and any nearby single target, school or ship position.

**Frequency difference.** Sv(higher frequency) minus Sv(lower frequency) on the MVBS grid, only where both have data (red = stronger at the higher frequency). It is useful for separating scatterer types, but only where both channels are above the noise. In the Baja recording 200 kHz is mostly noise below about 10 m, so the difference there is dominated by noise.

**Display controls.** The **Reset display** button returns the colour ranges to automatic, releases the lock, and resets the threshold, difference range, smoothing, depth and plot height. It leaves the cruise, recording, view and overlay checkboxes alone. *Hide Sv below* hides weak echoes (display only). Overlays toggle the seabed line, detected schools/layers (dashed white boxes), single targets (red circles) and layer boundaries.

**Other tabs.**

* **Regions:** hold Shift and drag on the MVBS echogram to draw boxes. The table gives cells, mean Sv, sa and NASC per frequency, plus the high-minus-low difference; type a class name per region and download the table (CSV) as annotations.
* **Integration:** the MVBS grid integrated into depth layers (5 to 50 m) for the whole recording or time bins, with coverage, mean Sv, sa, NASC and a count of single targets per layer.
* **Targets & Schools:** single-target table, TS histogram and TS-versus-depth plot, and the detected schools/layers.
* **Calibration & QC:** the calibration and environment values applied per channel (with a warning if built-in defaults were used) and data-quality checks (ping gaps, excluded pings, valid fraction, noise level).
* **Processing:** the processing steps in order, with the values used for this build.
* **Map:** the ship track is drawn when `nav.csv` exists.
* Every table has a CSV download button.

## Outputs (DuckDB tables)

| Table | One row per | Main columns |
|---|---|---|
| `echogram_mvbs` | file x channel x ping-time bin x depth bin | `frequency_khz`, `ping_time`, `depth_m`, `mvbs_db` |
| `echogram_seabed` | ping | `ping_time`, `seabed_depth_m`, `frequency_khz` (detection channel) |
| `echogram_files` | raw file | times, ping count, frequencies, waveform mode, sound speed, seabed range, `station_id`, `latitude`, `longitude`, SHA-256 |
| `echogram_stations` | survey in `stations.csv` | `station_id`, `latitude`, `longitude`, `position_confidence`, `processed` (a raw file matched) |
| `echogram_calibration` | file x channel | environment source, sound speed, absorption, gain, Sa correction, beam angle and widths, power, pulse length, impedances |
| `echogram_qc` | file x frequency | ping gaps, excluded pings, valid fraction above seabed, noise level |
| `echogram_targets` | detected single target | time, depth, TS (raw and beam-compensated), SNR, peak-over-background, along/athwart angle, echo length ratio |
| `echogram_schools` | detected school or layer | start/end time, depth range, thickness, mean and max Sv |
| `echogram_track` | ping (only with `nav.csv`) | latitude, longitude, distance, speed |
| `build_metadata` | build | config and raw-manifest SHA-256, Echopype version, git commit, row counts |

With `--write-parquet` the same tables are also written as Parquet under `processed/<cruise_id>/` (not written by default). Example query on the database:

    SELECT frequency_khz, depth_m, AVG(mvbs_db) FROM echogram_mvbs GROUP BY 1, 2;

### Zarr stores (full-resolution arrays)

The DuckDB tables hold the coarse MVBS grid and the detections. The full-resolution arrays are too large for a database, so they are written as Zarr stores (chunked arrays in a folder) under `processed/<cruise_id>/`, one per raw file:

| Store | Contents | Used for |
|---|---|---|
| `sv/<file_key>.zarr` | Cleaned, calibrated Sv at every ping and depth sample, all frequencies | The *Full resolution (cleaned Sv)* view, smoothing, automatic colour ranges and the Sv values in the details panel when you click the echogram |
| `raw_level/<file_key>.zarr` | Uncalibrated echo level at every ping and depth sample, all frequencies | The *Uncalibrated echo level* view |

The dashboard reads only the recording and frequency being shown, so these views stay fast. The paths are stored in the `echogram_files` table (`sv_zarr_path`, `raw_level_zarr_path`). Everything else (MVBS grid, difference view, seabed, targets, schools, integration, tables) comes from DuckDB and works without the Zarr stores. `--no-sv-store` skips the Sv stores; the output is much smaller, but the full-resolution view and smoothing are then unavailable.

To open a store in Python: `xr.open_zarr('processed/<cruise_id>/sv/<file_key>.zarr', consolidated=False)`.

## config.toml Schema

Section names are for readability only; every key must be unique.

### Environment (optional)
`TEMPERATURE_C`, `SALINITY_PSU`, `PRESSURE_DBAR` (set together; sound speed is derived with Mackenzie 1981) or `SOUND_SPEED_MS`. If absent, the values stored in the `.raw` file are used.

### Calibration (optional)
`ECS_FILE`: path to an EK80 calibration file from a sphere calibration. This code path is wired to Echopype's `ecs_file` option but has not been exercised on real data yet.

### Platform
`TRANSDUCER_DEPTH_M`: transducer depth below the surface, added to range to give depth.

### Processing
| Key | Default | Meaning |
|---|---|---|
| `NEAR_FIELD_EXCLUSION_M` | 3.0 | Samples closer than this are removed |
| `BACKGROUND_NOISE`, `BN_PING_NUM`, `BN_RANGE_M`, `SNR_THRESHOLD_DB` | true, 20, 5.0, 3.0 | Noise estimated per block of 20 pings x 5 m; cells under 3 dB SNR removed |
| `IMPULSE_NOISE`, `IMPULSE_DEPTH_BIN_M`, `IMPULSE_SIDE_PINGS`, `IMPULSE_THRESHOLD_DB` | true, 5.0, 2, 10.0 | Spike removal |
| `SEABED_CHANNEL_KHZ`, `SEABED_MIN_RANGE_M`, `SEABED_THRESHOLD_DB`, `SEABED_OFFSET_M`, `SEABED_MASK_MARGIN_M` | 38, 5.0, -40.0, 0.5, 1.0 | Seabed detector and mask |

### Excluded pings and transient noise
`EXCLUDE_TIME_RANGES_UTC` (list of `["start", "end"]` pairs) blanks bad pings. `TRANSIENT_NOISE` (default false), `TRANSIENT_DEPTH_BIN_M`, `TRANSIENT_SIDE_PINGS`, `TRANSIENT_THRESHOLD_DB`, `TRANSIENT_EXCLUDE_ABOVE_M` control the transient-noise mask (a lightweight version of the Ryan et al. 2015 method written for this project; Echopype's own was too slow and memory-hungry).

### Single targets (`[targets]`)
`TARGET_DETECTION` (true), `TARGET_TS_MIN_DB` (-60), `TARGET_SNR_MIN_DB` (10), `TARGET_MIN_LENGTH` / `TARGET_MAX_LENGTH` (0.8 / 1.8 pulse lengths), `TARGET_MAX_BEAM_COMP_DB` (6), `TARGET_MIN_DEPTH_M` (5), `TARGET_SEABED_MARGIN_M` (2), `TARGET_BACKGROUND_WINDOW_M` (5), `TARGET_PEAK_OVER_BACKGROUND_DB` (15; an echo must also stand this far above the local background, which keeps noise from being counted).

### Schools (`[schools]`)
`SCHOOL_DETECTION` (true), `SCHOOL_FREQUENCIES_KHZ` ([38]), `SCHOOL_THRESHOLD_DB` (-70), `SCHOOL_MIN_THICKNESS_M` (2), `SCHOOL_MIN_LENGTH_S` (5), `SCHOOL_SMOOTH_RANGE_M` (0.5), `SCHOOL_SMOOTH_PINGS` (5).

### Output
`MVBS_RANGE_BIN_M` (1.0) and `MVBS_PING_TIME_BIN` ("10s").

## Known limitations

* **Environment values.** Calibration depends on the water's temperature, salinity, pressure and sound speed, and these can come from CTD data. Put the values from a CTD cast taken close to the echosounder recording, in space and in time, in the cruise's `config.toml`; the farther the cast is from the recording, the less well it describes the water that was surveyed. Without an `[environment]` section the `.raw` file stores 10 C, 35 PSU and 100 m, the EK80 software defaults; the Calibration & QC tab warns when defaults were used.
* **Factory calibration only.** The stored gains are 18.0 dB (38 kHz) and 18.7 dB (200 kHz). No sphere-calibration file was supplied, so absolute Sv levels are uncertain. If you have an `.ecs` file, set `ECS_FILE`.
* **200 kHz noise.** In some recordings the 200 kHz channel is mostly noise below about 20 m: calibrated Sv rises steadily with range, the signature of a noise floor. The seabed is still visible. Treat 200 kHz values below the near surface with caution, and raise `SNR_THRESHOLD_DB` if you want stricter noise removal.
* **Positions are per survey, not per ping.** The `.raw` file contained no valid GPS data, so each file gets the single position of its survey from `stations.csv`. A moving survey would need the ship's navigation log joined by time.
* **Seabed detection.** A simple threshold detector. Check the red line on the echogram for each new cruise and adjust `SEABED_THRESHOLD_DB`.
* **Single targets are candidates.** Beam compensation needs split-beam angles, which Echopype provides for 38 kHz only here (200 kHz is single-beam), so 200 kHz TS is uncompensated. Without a net or camera, a detected target is not a confirmed fish. A recording can yield only a few 38 kHz candidates at the same depth, which looks more like a thin layer than fish.
* **Not implemented:** sphere calibration, per-ping impedance traces (impedance is stored once per file), video, fish tracking, multiple detection configurations, hand-edited seabed lines, real-time operation.
* **Sequential pinging.** If the channels transmit in turn (the Echoview tutorial data does), each channel's pings are paired in order and treated as simultaneous; the real time difference is about one ping interval and is ignored. A channel's extra unpaired ping at the end is dropped. Broadband (FM) files are calibrated in broadband mode with pulse compression.
* **Mixed CW and broadband files** are rejected with an error; split them first.

## Dashboard stops responding

The dashboard stops responding after a period of inactivity and needs to be refreshed (press F5 or reload the page). This is the browser's tab-sleeping setting, on Windows and on Mac: Edge puts an inactive tab to sleep (15 minutes on the Windows machine this was tested on; the time is configurable), and Chrome's Memory Saver can do the same. Either one cuts the tab's link to the dashboard server. In Edge, add the dashboard address to Settings > System and performance > Performance > "Always keep these sites active" (for example `localhost`), or turn off "Automatically put tabs to sleep". In Chrome, use Settings > Performance > "Always keep these sites active", or turn off Memory Saver.

## Contact

Joe Carlisle - joeacarlisle@gmail.com
