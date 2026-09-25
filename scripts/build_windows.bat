@echo off
rem Build of T-EVOLVE with MinGW-w64 gcc (MSYS2: pacman -S mingw-w64-ucrt-x86_64-gcc).
rem Usage: scripts\build_windows.bat      (run from the repository root)
setlocal
where gcc >nul 2>&1
if errorlevel 1 (
  echo gcc was not found. Install MSYS2 from https://www.msys2.org, then in the "MSYS2 UCRT64" shell run
  echo   pacman -S mingw-w64-ucrt-x86_64-gcc
  echo and add C:\msys64\ucrt64\bin to PATH.
  exit /b 1
)
if not exist build mkdir build
set COMMIT=unknown
for /f %%i in ('git rev-parse --short HEAD 2^>nul') do set COMMIT=%%i
set CF=-O3 -march=native -std=gnu11 -Wall -Wextra -Wno-unused-function -Iinclude -DGIT_COMMIT=\"%COMMIT%\" -DBUILD_FLAGS=\"gcc_-O3_-march=native_mingw\" -D__USE_MINGW_ANSI_STDIO=1
set CORE=src\keccak.c src\ring.c src\sample.c
set SRC=%CORE% src\params.c src\bdlop.c src\shamir.c src\codec.c src\ballot.c src\agg.c src\tally.c src\evolve.c
gcc %CF% -o build\test_core.exe tests\test_core.c %CORE% -lm || exit /b 1
gcc %CF% -o build\test_protocol.exe tests\test_protocol.c %SRC% -lm || exit /b 1
gcc %CF% -o build\exp_ballot_rej.exe experiments\exp_ballot_rej.c %SRC% -lm || exit /b 1
gcc %CF% -o build\bench_protocol.exe bench\bench_protocol.c %SRC% -lm || exit /b 1
gcc %CF% -o build\param_report.exe tools\param_report.c %SRC% -lm || exit /b 1
echo build ok (commit %COMMIT%)
endlocal
