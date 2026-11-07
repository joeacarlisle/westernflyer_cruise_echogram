@echo off
setlocal EnableDelayedExpansion
REM Rebuild every cruise under cruises into the echogram DuckDB database, one after another.
REM Usage:   echogram_build_all.bat [extra args forwarded to main.py, e.g. --range-bin 0.5]
REM Each cruise is a folder under cruises\ whose name is the cruise ID.
REM Add a new cruise by creating its folder; no edit to this script is needed.

set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%..\.."

if exist ".venv\Scripts\activate.bat" (
    call .venv\Scripts\activate.bat
) else (
    echo [WARNING] Virtual environment not found at .venv\Scripts\activate.bat
)

set "FAILED="
for /d %%C in (cruises\*) do (
    echo.
    echo ===== Building %%~nxC =====
    python main.py %%~nxC %*
    if errorlevel 1 set "FAILED=!FAILED! %%~nxC"
)

echo.
if defined FAILED (
    echo Build finished with failures:!FAILED!
    echo Check the newest files in logs\.
) else (
    echo All cruises built successfully.
)

popd
pause
