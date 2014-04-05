rem To run this file from CygWin prompt use: cmd /c BuildBindings.bat -u



rem set PATH=""
rem set LIBPATH=""
rem set INCLUDE=""
rem set LIB=""

rem c:
rem cd "c:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin"
rem call vcvars32.bat

rem p:

rem c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources C:\WPyRosetta\LinuxPyRosetta\rosetta.16384 %*

rem c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources P:\WinPyRosetta\rosetta.all %*



rem c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources P:\WinPyRosetta\rosetta %*



rem 64Bit version
rem use cmd /c BuildBindings.bat
rem %comspec% /k
call "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" amd64
c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources rosetta.window %*
