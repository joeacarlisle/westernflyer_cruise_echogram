#!/bin/bash
#
# Generic build script for any cruise.
#
# Usage:   ./echogram_build.sh <cruise_id>
# Example: ./echogram_build.sh 20250416_Baja
#
# Make executable once with:  chmod +x echogram_build.sh

set -u

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <cruise_id>"
    echo "Example: $0 20250416_Baja"
    read -p "Press any key to continue..."
    exit 2
fi

CRUISE_ID="$1"
shift  # any extra args are forwarded to main.py (e.g. --range-bin 0.5)

# Move to the project root (up two levels from where this script is located)
cd "$(dirname "$0")/../.."

# Activate the virtual environment
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
else
    echo "[WARNING] Virtual environment not found at .venv/bin/activate"
fi

# Note: build logs are now per-build/timestamped (logs/wf_echogram_build_<cruise>_<ts>.log),
# so we no longer delete the prior log.

echo "Starting Echogram build for cruise: ${CRUISE_ID}"
python main.py "${CRUISE_ID}" "$@"
status=$?

if [ "${status}" -eq 0 ]; then
    echo
    echo "Build finished successfully."
    echo "Logs are under logs/  (newest file is this run)."
else
    echo
    echo "Build failed (exit code ${status}). Check the newest file in logs/."
fi

read -p "Press any key to continue..."
exit ${status}
