@echo off
setlocal
wsl.exe -d Ubuntu --cd "%~dp0." -- bash ./run.sh
if errorlevel 1 pause
