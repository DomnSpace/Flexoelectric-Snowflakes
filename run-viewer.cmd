@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel% equ 0 (
  start "Snowflake Observatory Server" cmd /k "py -3 -m http.server 8088"
) else (
  where python >nul 2>nul
  if errorlevel 1 (
    echo Python 3 is required. Install Python 3 and retry.
    pause
    exit /b 1
  )
  start "Snowflake Observatory Server" cmd /k "python -m http.server 8088"
)
timeout /t 2 /nobreak >nul
start "" "http://localhost:8088/"
endlocal
