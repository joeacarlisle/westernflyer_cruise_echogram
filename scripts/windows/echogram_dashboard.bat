@echo off
:: Move to the project root (up two levels from where this script is located)
pushd "%~dp0..\.."

:: Check if the virtual environment exists in the root
if exist ".venv\Scripts\python.exe" (
    echo [INFO] Found virtual environment. Launching...
    set PYTHON_EXE=.venv\Scripts\python.exe
) else (
    echo [INFO] Virtual environment not found. Attempting to use global Python...
    set PYTHON_EXE=python
)

:: Run the Panel dashboard
:: We use '-m panel' to ensure it uses the version installed in the environment
%PYTHON_EXE% -m panel serve echogram_holoviews.py --port 5008 --show

:: Keep the window open if the script crashes
if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] The dashboard failed to launch. Please check the logs above.
    pause
)

:: Return to original directory
popd
