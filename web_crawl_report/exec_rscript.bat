@ECHO OFF

setlocal

set "psCommand="(new-object -COM 'Shell.Application').BrowseForFolder(0,'Please choose a folder.',0,0).self.path""

for /f "usebackq delims=" %%f in (`powershell %psCommand%`) do set "file_path=%%f"

IF NOT DEFINED file_path (
	echo ERROR: You must choose a directory^!
	pause
	exit /b 1 
)

setlocal enabledelayedexpansion

Rscript %~dp0\report_qc.R !file_path! %~dp0\

endlocal

echo FINISHED.
PAUSE


@ECHO OFF

setlocal

set "psCommand="(new-object -COM 'Shell.Application').BrowseForFolder(0,'Please choose a folder.',0,0).self.path""

for /f "usebackq delims=" %%f in (`powershell %psCommand%`) do set "file_path=%%f"

IF NOT DEFINED file_path (
	echo ERROR: You must choose a directory^!
	pause
	exit /b 1 
)

setlocal enabledelayedexpansion

Rscript %~dp0\gnomad_af_report.R !file_path! %~dp0\

endlocal

echo FINISHED.
PAUSE