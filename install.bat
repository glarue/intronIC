@echo off
REM
REM Simple installation script for intronIC (Windows)
REM Usage: install.bat [--dev]
REM

echo ======================================
echo intronIC Installation Script
echo ======================================
echo.

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: python not found. Please install Python 3.8 or later.
    echo Download from: https://www.python.org/downloads/
    pause
    exit /b 1
)

REM Display Python version
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo Found Python %PYTHON_VERSION%

REM Check if pip is available
python -m pip --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: pip not found. Please install pip.
    pause
    exit /b 1
)

REM Determine installation mode
if "%1"=="--dev" (
    echo Installing in DEVELOPMENT mode...
    python -m pip install -e ".[dev]"
) else (
    echo Installing in USER mode...
    echo (Use 'install.bat --dev' for development mode)
    python -m pip install .
)

echo.
echo ======================================
echo Installation Complete!
echo ======================================
echo.
echo Verify installation:
echo   intronIC --version
echo   intronIC --help
echo.
echo Run a test:
echo   intronIC -g genome.fa -a annotation.gff3 -n species_name
echo.
echo For more information, see:
echo   - INSTALL.md for detailed installation options
echo   - README.md for usage guide
echo   - intronIC --help for all command-line options
echo.
pause
