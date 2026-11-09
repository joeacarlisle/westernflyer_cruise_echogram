# Project Components and Architecture

This document describes the functional components that comprise the Western Flyer EK80 echogram processing pipeline and provides guidelines for AI agents interacting with this repository.

## Functional Components

The codebase is organized into modular components, each responsible for a distinct stage of the data lifecycle.

### 1. Ingestion Component (`ek80_ingestion.py`)
**Role:** Raw Data Interface & Calibration
- **Responsibilities:**
    - Reads Simrad EK80 `.raw` files with Echopype and, only when asked (`--write-netcdf`), writes a lossless SONAR-netCDF4 archive (`.nc`).
    - Loads cruise-specific parameters from `config.toml` and survey positions from `stations.csv`.
    - Chooses the calibration mode from what was actually transmitted (CW vs broadband; complex vs power).
    - Applies optional environment overrides (temperature, salinity, pressure, sound speed) and an optional `.ecs` calibration file.
    - Calibrates to Sv (dB re 1 m-1) and adds a depth axis.
- **Key Interfaces:** `ingest_ek80_raw()`, `calibrate_sv()`, `load_config()`, `load_stations()`

### 2. Processing Component (`echogram_processing.py`)
**Role:** Echogram Cleaning & Gridding
- **Responsibilities:**
    - Near-field exclusion, background-noise removal (De Robertis & Higginbottom 2007), impulse-noise masking (Ryan et al. 2015).
    - Seabed detection and masking below the seabed.
    - Gridding to MVBS (mean volume backscattering strength) cells.
    - Flattening MVBS and seabed results to tables for Parquet/DuckDB.
    - All steps run **per channel**, because the two EK80 channels are sampled at different rates.
    - Excluded-ping and transient-noise handling, QC summary (`qc_summary()`), school detection (`find_schools()`).
- **Key Interfaces:** `process_sv()`, `compute_mvbs()`, `validate_config()`

### 2b. Analysis helpers (`echogram_analysis.py`, `echogram_targets.py`)
**Role:** Pure computation shared by the build and the dashboard (no I/O).
- `echogram_analysis.py`: block averaging, `integrate_layers()` (sa, NASC), `region_stats()`, `frequency_difference()`, `detect_schools()`, `load_nav()`/`interpolate_track()`, `transient_noise_mask()`.
- `echogram_targets.py`: `detect_single_targets()` (TS, beam compensation, peak-over-background test).

### 3. Orchestrator Component (`main.py`)
**Role:** Workflow Manager & Persistence
- **Responsibilities:**
    - Orchestrates the end-to-end pipeline for a specific cruise.
    - Manages CLI arguments, logging, and directory resolution.
    - Joins survey positions onto each file by `station_id` (`attach_positions()`).
    - Writes cleaned-Sv Zarr stores and atomic, idempotent DuckDB commits. The `.nc` archive and the Parquet copies are written only with `--write-netcdf` and `--write-parquet`. A table whose columns changed is recreated, so rebuild other cruises afterwards.
    - Records build provenance (git commit hash, `config.toml` SHA-256, raw-file SHA-256 manifest, Echopype version, timestamps).
- **Key Interface:** `main()` entry point and `process_one_file()`.

### 4. Analysis Component (`echogram_holoviews.py`)
**Role:** Interactive Visualization
- **Responsibilities:**
    - Serves a Panel dashboard (tabs: Echogram with cursor details and overlays, Regions, Integration, Targets & Schools, Seabed, Map, Calibration & QC, Processing, MVBS table, Cruise summary; CSV download on tables): full-resolution echogram (Datashader-rasterized), MVBS grid view, an uncalibrated echo-level view, one panel per detected frequency with per-frequency automatic colour ranges and an optional lock, an adjustable plot height, optional block-average smoothing of the full-resolution view, seabed depth, a map of survey positions, MVBS table, cruise summary and build provenance.
    - Reads small tables from DuckDB and the full-resolution Sv lazily from Zarr, so large recordings never go to the browser sample by sample.
- **Key Interface:** `dashboard.servable()`

### 5. Infrastructure Component (Scripts)
**Role:** Environment-Specific Automation
- **Linux/macOS (`scripts/linux_mac/`):** Bash scripts for building the database and launching the dashboard.
- **Windows (`scripts/windows/`):** Batch files for equivalent operations in a CMD environment.

---

## AI Agent Interaction Guidelines

When an AI agent interacts with this repository, it should adhere to the following operational patterns:

### Adding a New Cruise
1. **Directory Structure:** Create `cruises/<cruise_id>/`.
2. **Raw Data:** Place `.raw` files in `cruises/<cruise_id>/raw/`.
3. **Configuration:** A cruise with no `config.toml` gets one written automatically at build time with default settings and `CONFIG_CHECKED = false` (or run `python main.py <cruise_id> --init-config`; or copy one from an existing cruise). While `CONFIG_CHECKED = false` the build logs REVIEW REQUIRED and the dashboard shows a red banner on the Calibration & QC and Processing tabs; the user sets it to true after checking the seabed line, noise removal and water values. Keys missing from a cruise file fall back to `config.defaults.toml`. Never put measured environment values in the defaults file.
4. **Positions:** If `stations.csv` is missing the build writes one from the GPS fixes inside the raw files (and a `nav.csv` if the ship moved); check it. Otherwise add `stations.csv` with a `station_id` (start of the raw file name, before the first `-`) and decimal-degree `lat` and `lon` (negative south/west) per survey. These three column names match the CTD project's `stations.csv`.

### Modifying the Pipeline
- **Calibration or ingestion changes:** `ek80_ingestion.py`. Use Echopype's public API (`open_raw`, `compute_Sv`) so results stay comparable with other Echopype users.
- **Cleaning, seabed or gridding changes:** `echogram_processing.py`. Keep every step per channel.
- **Database Schema:** If adding columns, update the table-building functions in `echogram_processing.py`, the commit logic in `main.py`, and the queries in `echogram_holoviews.py`.
- Do not modify the sibling CTD project (`westernflyer_cruise_ctd`); this repository is independent of it.

### Execution and Verification
- **Run the Build:** Use `python main.py <cruise_id>` to process data.
- **Logs:** Monitor `logs/wf_echogram_build_<cruise_id>_<timestamp>.log` for detailed processing diagnostics.
- **Database:** Verify output in `processed/wf_echogram.duckdb`. The `build_metadata` and `echogram_files` tables should be checked to confirm successful persistence and provenance.

### Dependencies
- This project requires **Python 3.13 or newer** (Echopype 0.12).
- Other dependencies are listed in `pyproject.toml`.
