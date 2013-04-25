rem set PATH=""
rem set LIBPATH=""
rem set INCLUDE=""
rem set LIB=""

rem c:
rem cd "c:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin"
rem call vcvars32.bat

rem p:

rem c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources C:\WPyRosetta\LinuxPyRosetta\rosetta.16384 %*

c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources P:\WinPyRosetta\rosetta %*
rem c:\Python27\python.exe BuildBindings.py --use-pre-generated-sources P:\WinPyRosetta\rosetta.all %*
