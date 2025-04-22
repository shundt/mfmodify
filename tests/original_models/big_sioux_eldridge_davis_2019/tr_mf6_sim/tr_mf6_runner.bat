@echo off
:: Batch file used to run the MF6 version of the Sioux Falls Big Sioux aquifer GW Model
:: Kyle W Davis (kyledavis@usgs.gov) November 13, 2017
:: updated November 14, 2017 - kyledavis


:: Setting variables
if not defined wkdir (goto :PATHS) else (goto :RUN)


:PATHS
cd ..\..\
set wkdir=%cd%

set moddir6=tr_mf6_sim
set modbase6=sfbsa_tr6
set namfile6=%modbase6%.nam
set mf=MODFLOW 6
set mf_ver=6.0.3
set mf_base=mf6
set mf_exe=%mf_base%.exe
set modnam=Sioux Falls Big Sioux Aquifer
set zbud6=zbud6.exe
set python=%wkdir%\ancillary\SFBSA_py36x64\python.exe
goto :RUN


:RUN
:: Running the model
if not defined call_runner (goto :MODINFO) else (goto :OUTCHEK)


:MODINFO
echo.
echo %modnam% GW model using using %mf% v%mf_ver%
DATE /t
TIME /t
echo.
echo Master model directory:
echo %wkdir%
cd %wkdir%\model\%moddir6%\
echo.
echo Deleting old results if they exist
if exist %wkdir%\output\output.%moddir6% RD %wkdir%\output\output.%moddir6% /S /Q
echo.
goto :OUTCHEK


:OUTCHEK
if not exist %wkdir%\output\output.%moddir6%\ (goto :MKDIRS) else (goto :MODRUN)


:MKDIRS
echo Making output directories if they do not exist
echo.
if not exist %wkdir%\output\ mkdir %wkdir%\output\
if not exist %wkdir%\output\output.%moddir6%\ mkdir %wkdir%\output\output.%moddir6%\
if not exist %wkdir%\output\output.%moddir6%\obs\ mkdir %wkdir%\output\output.%moddir6%\obs\
goto :MODRUN


:MODRUN
:: NOTE...modflow6 v6.0.3 does not need the name file referenced after the .exe (previous versions do need it)
set call_zbud=True
%wkdir%\bin\mf6.exe
:: copy outputs
robocopy .\ %wkdir%\output\output.%moddir6% mfsim.lst *.grb /XO /MOV > NUL
:: copy grid discretization file
robocopy %wkdir%\model\external_files\ %wkdir%\output\output.%moddir6% *.grb /XO /MOV > nul
:: run zonebudget
cd %wkdir%\ancillary\zonebudget\tr_mf6_sim\
call | .\tr_zbud6_process.bat

if not defined call_runner (goto :PRINT) else (goto :EOF)


:PRINT
echo Finished!
echo.
echo Results in:
echo %wkdir%\output\output.%moddir6%\
echo.
pause
goto :EOF
