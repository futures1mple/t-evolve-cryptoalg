@echo off
rem Scaling measurements on Windows: ballots for n = 4..15 and L = 1..10, aggregation for
rem N_V = 10^4 (3 runs), 5*10^4 and 10^5 (2 runs), all in the low-memory sequential mode.
rem
rem   scripts\run_windows_scaling.bat          everything (about 3.5-4 hours)
rem   scripts\run_windows_scaling.bat ballots  only the ballots (about 5 minutes)
rem   scripts\run_windows_scaling.bat agg      only aggregation 10^4 and 5*10^4 (about 1 hour)
rem   scripts\run_windows_scaling.bat big      only aggregation 10^5 (about 2-2.5 hours, needs about 11 GB of free memory)
rem
rem Before running: plug in the power, choose the best-performance power mode, close other
rem programs, disable sleep. Results go to results\windows\scaling_<date>.
setlocal
set PART=%1
if "%PART%"=="" set PART=all
call scripts\build_windows.bat || exit /b 1
for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmm"') do set DT=%%i
set OUT=results\windows\scaling_%DT%
if not exist %OUT% mkdir %OUT%
echo writing to %OUT%
powershell -NoProfile -Command "Get-CimInstance Win32_Processor | Format-List Name,NumberOfCores,MaxClockSpeed; 'TotalPhysicalMemoryGB: ' + [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1); 'FreeMemoryGB: ' + [math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB,1)" > %OUT%\machine.txt 2>&1
powercfg /getactivescheme >> %OUT%\machine.txt 2>&1
(gcc --version & ver) >> %OUT%\machine.txt 2>&1
echo part: %PART% >> %OUT%\machine.txt
echo [check] core tests
build\test_core.exe > %OUT%\test_core.txt || (echo core tests FAILED & exit /b 1)

if "%PART%"=="all" goto ballots
if "%PART%"=="ballots" goto ballots
goto afterballots
:ballots
echo [ballots 1/2] yes/no ballots, n = 4, 5, 7, 10, 15 (q = 2^42 - 143)
build\bench_protocol.exe -reps 30 -ballot 4,3,1,-1 -ballot 5,3,1,-1 -ballot 7,4,1,-1 -ballot 10,6,1,-1 -ballot 15,8,1,-1 -out %OUT%\ballots.csv > %OUT%\ballots_stdout.txt || (echo ballots FAILED & exit /b 1)
echo [ballots 2/2] one of L candidates, L = 2, 5 (q = 2^42 - 143) and L = 10 (q = 2^43 - 175)
build\bench_protocol.exe -reps 30 -ballot 4,3,2,1 -ballot 4,3,5,1 -ballot 10,6,2,1 -out %OUT%\ballots.csv >> %OUT%\ballots_stdout.txt || (echo ballots FAILED & exit /b 1)
build\bench_protocol.exe -reps 30 -q 8796093022033 -ballot 4,3,10,1 -out %OUT%\ballots.csv >> %OUT%\ballots_stdout.txt || (echo ballots FAILED & exit /b 1)
if "%PART%"=="ballots" goto done
:afterballots

if "%PART%"=="all" goto agg
if "%PART%"=="agg" goto agg
goto afteragg
:agg
echo [agg 1/2] aggregation, 10^4 voters, 3 runs (about 20 minutes)
for %%s in (1 2 3) do build\bench_protocol.exe -only-agg -seq -nv 10000 -seed %%s -out %OUT%\agg.csv >> %OUT%\agg_stdout.txt
echo [agg 2/2] aggregation, 5*10^4 voters (about 40 minutes)
build\bench_protocol.exe -only-agg -seq -nv 50000 -seed 1 -out %OUT%\agg.csv >> %OUT%\agg_stdout.txt
if "%PART%"=="agg" goto done
:afteragg

if "%PART%"=="all" goto big
if "%PART%"=="big" goto big
goto done
:big
echo [big] aggregation, 10^5 voters, 2 runs (the longest step; needs about 11 GB of free memory)
for %%s in (1 2) do build\bench_protocol.exe -only-agg -seq -nv 100000 -seed %%s -out %OUT%\agg_big.csv >> %OUT%\agg_big_stdout.txt

:done
echo done. Please send the folder %OUT%
endlocal
