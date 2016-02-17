@call :GetRDir
@if "%R%"=="" goto error_no_RDIR

@call :GetRToolsDir
@if "%RTOOLS%"=="" goto error_no_RTOOLSDIR

@if not "%*"=="" SET MinGW=%*
@if "%MinGW%"=="" call :GetMinGWDir
@REM if "%MinGW%"=="" goto error_no_MinGWDIR

@rem PATH
@rem ----
@if exist "%VSINSTALLDIR%Team Tools\Performance Tools" @set PATH=%VSINSTALLDIR%Team Tools\Performance Tools;%PATH%
@if exist "%ProgramFiles%\HTML Help Workshop" set PATH=%ProgramFiles%\HTML Help Workshop;%PATH%
@if exist "%ProgramFiles(x86)%\HTML Help Workshop" set PATH=%ProgramFiles(x86)%\HTML Help Workshop;%PATH%
@if exist "%VCINSTALLDIR%VCPackages" set PATH=%VCINSTALLDIR%VCPackages;%PATH%
@if exist "%FrameworkDir%%Framework40Version%" set PATH=%FrameworkDir%%Framework40Version%;%PATH%
@if exist "%FrameworkDir%%FrameworkVersion%" set PATH=%FrameworkDir%%FrameworkVersion%;%PATH%
@if exist "%VSINSTALLDIR%Common7\Tools" set PATH=%VSINSTALLDIR%Common7\Tools;%PATH%
@if exist "%VCINSTALLDIR%BIN" set PATH=%VCINSTALLDIR%BIN;%PATH%
@set PATH=%R%\bin;%RTOOLS%\gcc-4.6.3\bin;%PATH%
@REM set PATH=%MinGW%\bin;%PATH%
@set PATH=%RTOOLS%\bin;%PATH%
@set PATH=%LOCALAPPDATA%\Pandoc;%PATH%

@goto end

@REM -----------------------------------------------------------------------
:GetRDir
@set R=
@call :GetRDirHelper32 HKLM > nul 2>&1
@if errorlevel 1 call :GetRDirHelper32 HKCU > nul 2>&1
@if errorlevel 1 call :GetRDirHelper64 HKLM > nul 2>&1
@if errorlevel 1 call :GetRDirHelper64 HKCU > nul 2>&1
@if errorlevel 1 call :GetRDirHelper
@if errorlevel 1 exit /B 1
@exit /B 0

:GetRDirHelper32
@for /F "tokens=1,2*" %%i in ('reg query "%1\SOFTWARE\R-core\R" /v "InstallPath"') DO (
	@if "%%i"=="InstallPath" (
		@SET R=%%k
	)
)
@if "%R%"=="" exit /B 1
@exit /B 0

:GetRDirHelper64
@for /F "tokens=1,2*" %%i in ('reg query "%1\SOFTWARE\Wow6432Node\R-core\R" /v "InstallPath"') DO (
	@if "%%i"=="InstallPath" (
		@SET R=%%k
	)
)
@if "%R%"=="" exit /B 1
@exit /B 0

:GetRDirHelper
@if exist "%ProgramFiles%\R" @SET R_ROOT=%ProgramFiles%\R
@for /f "delims=" %%a in (
    'dir /b /od /ad "%R_ROOT%" 2^>NUL'
) do @SET R=%R_ROOT%\%%a
@if "%R%"=="" exit /B 1
@exit /B 0

@REM -----------------------------------------------------------------------
:GetRToolsDir
@set RTOOLS=
@call :GetRToolsDirHelper32 HKLM > nul 2>&1
@if errorlevel 1 call :GetRToolsDirHelper32 HKCU > nul 2>&1
@if errorlevel 1 call :GetRToolsDirHelper64 HKLM > nul 2>&1
@if errorlevel 1 call :GetRToolsDirHelper64 HKCU > nul 2>&1
@if errorlevel 1 call :GetRToolsDirHelper
@if errorlevel 1 exit /B 1
@exit /B 0

:GetRToolsDirHelper32
@for /F "tokens=1,2*" %%i in ('reg query "%1\SOFTWARE\R-core\Rtools" /v "InstallPath"') DO (
	@if "%%i"=="InstallPath" (
		@SET RTOOLS=%%k
	)
)
@if "%RTOOLS%"=="" exit /B 1
@exit /B 0

:GetRToolsDirHelper64
@for /F "tokens=1,2*" %%i in ('reg query "%1\SOFTWARE\Wow6432Node\R-core\Rtools" /v "InstallPath"') DO (
	@if "%%i"=="InstallPath" (
		@SET RTOOLS=%%k
	)
)
@if "%RTOOLS%"=="" exit /B 1
@exit /B 0

:GetRToolsDirHelper
@if exist "C:\Rtools\bin" @SET RTOOLS=c:\Rtools
@if "%RTOOLS%"=="" exit /B 1
@exit /B 0

@REM -----------------------------------------------------------------------
:GetMinGWDir
@set MinGW=
@if exist "C:\MinGW\bin\gcc.exe" @SET MinGW=C:\MinGW
@if "%MinGW%"=="" exit /B 1
@exit /B 0

@REM -----------------------------------------------------------------------
:error_no_RDIR
@echo ERROR: Cannot determine the location of the R folder.
@goto end

:error_no_RTOOLSDIR
@echo ERROR: Cannot determine the location of the RTools folder.
@goto end

:error_no_MinGWDIR
@echo ERROR: Cannot determine the location of the MinGW folder.
@goto end

:end
