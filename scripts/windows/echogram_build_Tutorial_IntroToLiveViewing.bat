@echo off
REM Convenience wrapper for the Tutorial_IntroToLiveViewing cruise (Echoview tutorial data, 120 and 200 kHz).
REM Other cruises: echogram_build.bat ^<cruise_id^>
call "%~dp0echogram_build.bat" Tutorial_IntroToLiveViewing %*
