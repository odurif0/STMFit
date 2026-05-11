@echo off
setlocal enabledelayedexpansion
title Multi-Gaussian Fit GUI

set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

set "JULIA="
for /f "delims=" %%i in ('where julia 2^>nul') do (
    if exist "%%i" ( set "JULIA=%%i" & goto :julia_found )
)
if exist "%LOCALAPPDATA%\Microsoft\WindowsApps\julia.exe" set "JULIA=%LOCALAPPDATA%\Microsoft\WindowsApps\julia.exe" & goto :julia_found
if exist "%PROGRAMFILES%\Julia-1.10\bin\julia.exe"   set "JULIA=%PROGRAMFILES%\Julia-1.10\bin\julia.exe"   & goto :julia_found
if exist "%PROGRAMFILES%\Julia-1.11\bin\julia.exe"   set "JULIA=%PROGRAMFILES%\Julia-1.11\bin\julia.exe"   & goto :julia_found

echo Julia is not installed. Installing via winget...
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
if errorlevel 1 (
    echo ERROR: Please install Julia manually from https://julialang.org/downloads/
    pause
    exit /b 1
)

:julia_found
echo Julia: %JULIA%
echo Installing dependencies...

"%JULIA%" -e "using Pkg; Pkg.develop(path=\"%USERPROFILE%\\Git\\GaussianFit1D.jl\"); Pkg.instantiate()" 2>nul
"%JULIA%" --project=. -e "using Pkg; Pkg.develop(path=\"%USERPROFILE%\\Git\\GaussianFit1D.jl\"); Pkg.instantiate()"

echo Starting GUI...
"%JULIA%" --project=. app.jl

pause
