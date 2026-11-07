@echo off
REM Generic build script for any cruise.
REM
REM Usage:   echogram_build.bat <cruise_id> [extra args forwarded to main.py]
REM Example: echogram_build.bat 20250416_Baja
REM Example: echogram_build.bat 20250416_Baja --range-bin 0.5

REM Capture script directory BEFORE any shift calls.
REM %~dp0 is derived from %0, which shifts along with %1..%9 on each
REM shift call. Storing it now keeps the correct path regardless of how
REM many arguments are consumed below.
set "SCRIPT_DIR=%~dp0"

if "%~1"=="" (
    echo Usage: %~nx0 ^<cruise_id^> [extra args]
    echo Example: %~nx0 20250416_Baja
    pause
    exit /b 2
)

set "CRUISE_ID=%~1"
shift

REM Collect any remaining args into EXTRA_ARGS (since %* does not update after shift)
set "EXTRA_ARGS="
:collect_args
if "%~1"=="" goto done_collect
set "EXTRA_ARGS=%EXTRA_ARGS% %~1"
shift
goto collect_args
:done_collect

REM Move to the project root (up two levels from scripts\windows\)
pushd "%SCRIPT_DIR%..\.."

REM Activate the virtual environment
if exist ".venv\Scripts\activate.bat" (
    call .venv\Scripts\activate.bat
) else (
    echo [WARNING] Virtual environment not found at .venv\Scripts\activate.bat
)

echo Starting Echogram build for cruise: %CRUISE_ID%
python main.py %CRUISE_ID%%EXTRA_ARGS%

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Build finished successfully.
    echo Logs are under logs\  (newest file is this run^).
) else (
    echo.
    echo Build failed ^(exit code %ERRORLEVEL%^). Check the newest file in logs\.
)

popd
pause
