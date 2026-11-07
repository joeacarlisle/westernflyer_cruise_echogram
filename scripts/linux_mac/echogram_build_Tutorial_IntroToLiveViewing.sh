#!/bin/bash
# Convenience wrapper for the Tutorial_IntroToLiveViewing cruise (Echoview tutorial data, 120 and 200 kHz).
# Other cruises: ./echogram_build.sh <cruise_id>
exec "$(dirname "$0")/echogram_build.sh" Tutorial_IntroToLiveViewing "$@"
