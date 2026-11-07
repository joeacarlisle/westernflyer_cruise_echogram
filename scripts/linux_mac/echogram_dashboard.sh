#!/bin/bash
#
# Make executable once with:  chmod +x echogram_dashboard.sh

# Move to the project root (up two levels from where this script is located)
cd "$(dirname "$0")/../.."

# Check if the virtual environment exists in the root (Unix paths use forward slashes)
if [ -f ".venv/bin/python" ]; then
    echo "[INFO] Found virtual environment. Launching..."
    PYTHON_EXE=".venv/bin/python"
else
    echo "[INFO] Virtual environment not found. Attempting to use global Python..."
    PYTHON_EXE="python"
fi

# Run the Panel dashboard
# We use '-m panel' to ensure it uses the version installed in the environment
# REMOVE '--show' if running on headless Linux servers as it attempts to open a browser window
"$PYTHON_EXE" -m panel serve echogram_holoviews.py --port 5008 --show

# Keep the window open if the script crashes
if [ $? -ne 0 ]; then
    echo ""
    echo "[ERROR] The dashboard failed to launch. Please check the logs above."
    read -p "Press any key to continue..."
fi
