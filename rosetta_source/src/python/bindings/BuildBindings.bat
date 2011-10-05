set PATH=""
set LIBPATH=""
set INCLUDE=""
set LIB=""

c:
cd "c:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin"
call vcvars32.bat

p:
c:\Python27\python.exe BuildBindings.py


