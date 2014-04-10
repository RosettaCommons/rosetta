@echo off

rem To run this file from CygWin prompt use: cmd /c BuildBindings.bat -u

rem c:
rem cd "c:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin"
rem call vcvars32.bat

rem c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources P:\WinPyRosetta\rosetta %*

rem 64Bit version
rem use cmd /c BuildBindings.bat
rem %comspec% /k
call "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" amd64
c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources rosetta.windows %*
