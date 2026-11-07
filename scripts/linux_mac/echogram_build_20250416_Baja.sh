#!/bin/bash
# Convenience wrapper for the 20250416_Baja cruise (Baja California 2025).
# Other cruises: ./echogram_build.sh <cruise_id>
exec "$(dirname "$0")/echogram_build.sh" 20250416_Baja "$@"
