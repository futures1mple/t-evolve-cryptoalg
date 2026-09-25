@echo off
rem Complete measurement run on Windows. Takes about one hour on a desktop CPU.
rem Before running: plug in the power, choose the "High performance" / "Best performance" power
rem mode, close other programs. Results go to results\windows\<date>.
setlocal
call scripts\build_windows.bat || exit /b 1
for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmm"') do set DT=%%i
set OUT=results\windows\%DT%
if not exist %OUT% mkdir %OUT%
echo writing to %OUT%
powershell -NoProfile -Command "Get-CimInstance Win32_Processor | Format-List Name,NumberOfCores,MaxClockSpeed; Get-CimInstance Win32_PhysicalMemory | Format-List Capacity,Speed; (Get-CimInstance Win32_PowerPlan -Namespace root\cimv2\power -Filter 'IsActive=True').ElementName" > %OUT%\machine.txt 2>&1
(gcc --version & ver) >> %OUT%\machine.txt 2>&1
echo [1/5] tests
build\test_core.exe > %OUT%\test_core.txt || (echo core tests FAILED & exit /b 1)
build\test_protocol.exe -q > %OUT%\test_protocol.txt || (echo protocol tests FAILED & exit /b 1)
echo [2/5] ballot rejection experiment (1000 ballots per configuration)
build\exp_ballot_rej.exe 1000 20260924 %OUT% > %OUT%\ballot_rej_stdout.txt
echo [3/5] micro and ballot benchmarks
build\bench_protocol.exe -reps 50 -skip-agg -out %OUT%\bench.csv > %OUT%\bench_stdout.txt
echo [4/5] aggregation, 1000 voters
build\bench_protocol.exe -only-agg -nv 1000 -out %OUT%\bench.csv >> %OUT%\bench_stdout.txt
echo [5/5] aggregation, 10000 voters (the longest step)
build\bench_protocol.exe -only-agg -nv 10000 -out %OUT%\bench.csv >> %OUT%\bench_stdout.txt
echo done. Please send the folder %OUT%
endlocal
