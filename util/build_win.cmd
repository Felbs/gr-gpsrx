@echo off
REM Build the C++ blocks on Windows against radioconda with MSVC Build Tools + the Ninja and
REM CMake that ship with Visual Studio. Run from the repo root:  util\build_win.cmd [install]
REM Assumes radioconda at %USERPROFILE%\radioconda and VS 2022 Build Tools.
setlocal
set VSROOT=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools
set CONDA=%USERPROFILE%\radioconda
call "%VSROOT%\VC\Auxiliary\Build\vcvars64.bat" >nul
set PATH=%VSROOT%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin;%VSROOT%\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja;%CONDA%\Library\bin;%CONDA%;%PATH%
if not exist build\cpp mkdir build\cpp
cd build\cpp
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=%CONDA:\=/%/Library -DPYTHON_EXECUTABLE=%CONDA:\=/%/python.exe -DCMAKE_INSTALL_PREFIX=%CONDA:\=/%/Library -DENABLE_DOXYGEN=OFF ../.. || exit /b 1
ninja || exit /b 1
if "%1"=="install" ninja install
