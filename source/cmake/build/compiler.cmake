## Note: This build uses the CMAKE default compiler settings for your computer.
## Changing the "COMPILER" setting here will *not* change which compiler is used;
## it only changes the extension of the executable (e.g. .linuxgccdebug).
##
## To change the default compiler for CMAKE, set the $CC and $CXX shell environment variables.
## For example, in bash, "export CXX=clang++; export CC=clang;"
## Or "export CXX='distcc g++'; export CC='distcc gcc';"
##
## Alternatively, you can pass the -DCMAKE_CXX_COMPILER="<compiler name>" flag
## when you invoke cmake.
##
## IMPORTANT: These setting are only obeyed on the first invocation on a clean build.
## An incremental build will not change the compiler even if the $CC/$CXX settings change.
## Because of this - and unlike the scons system - the cmake system can't handle
## multiple compilers on the same target. If you need multiple simultaneous compiler
## builds, use additional build directories.


## Try to set the executable extension appropriately for people with a default compiler of clang.
## Assume that if you're not obviously running clang, you're running gcc.
## There's additional CMAKE flags that can tease other cases apart, but
## this /should/ work for most of Rosetta developers.

IF( ${CMAKE_CXX_COMPILER_ID} MATCHES ".*[Cc][Ll][Aa][Nn][Gg].*" )
    SET(COMPILER clang)
ELSE()
    SET(COMPILER gcc)
ENDIF()

MESSAGE( ">> CMAKE identifies C++ compiler as '${CMAKE_CXX_COMPILER_ID}', interpreting this as '${COMPILER}'" )
MESSAGE( ">> To change, set CXX and CC environment variables (or pass -DCMAKE_CXX_COMPILER) and do a clean rebuild.")
MESSAGE( ">>  current settings: CXX='$ENV{CXX}' CC='$ENV{CC}'" )


MESSAGE(STATUS "Compiler: ${COMPILER}")



