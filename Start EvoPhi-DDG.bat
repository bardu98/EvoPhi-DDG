@echo off
setlocal
wsl.exe --status >nul 2>&1
if errorlevel 1 (
  echo This version requires WSL2 with a Linux distribution installed.
  echo See gui/README.md for Windows setup.
  pause
  exit /b 1
)
for /f "delims=" %%p in ('wsl.exe wslpath -u "%~dp0."') do set "EVOPHI_REPO=%%p"
if not defined EVOPHI_REPO (
  echo Could not locate the project in WSL. See gui/README.md.
  pause
  exit /b 1
)
wsl.exe bash "%EVOPHI_REPO%/gui/start.sh"
if errorlevel 1 pause
