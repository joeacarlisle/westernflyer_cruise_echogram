#!/bin/bash
#
# Rebuild every cruise under cruises/ into the DuckDB database, one after another.
#
# Usage:   ./echogram_build_all.sh [extra args forwarded to main.py, e.g. --range-bin 0.5]
# Each cruise is a folder under cruises/ whose name is the cruise ID.
# Add a new cruise by creating its folder; no edit to this script is needed.
#
# Make executable once with:  chmod +x echogram_build_all.sh

cd "$(dirname "$0")/../.."

if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
else
    echo "[WARNING] Virtual environment not found at .venv/bin/activate"
fi

failed=()
for dir in cruises/*/; do
    cruise="$(basename "$dir")"
    echo
    echo "===== Building ${cruise} ====="
    python main.py "${cruise}" "$@" || failed+=("${cruise}")
done

echo
if [ "${#failed[@]}" -gt 0 ]; then
    echo "Build finished with failures: ${failed[*]}"
    echo "Check the newest files in logs/."
    read -p "Press any key to continue..."
    exit 1
fi
echo "All cruises built successfully."
read -p "Press any key to continue..."
