CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

##############################################################################
# This file only needs to specify which mode (i.e., release, debug),
# extras (i.e., graphics), and compiler (i.e., gcc, clang) should be used. The
# appropriate settings for each mode, extras, and compiler should be
# implemented in ../build/build.settings.cmake.
##############################################################################

############ Build Type, Mode, Link Rules, Compiler, Extras
SET(MODE release)
SET(COMPILER gcc)
SET(LINK_RULES shared)

############ Build Rules
INCLUDE(../build/compiler.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)
# INCLUDE(../build/build.settings.cmake)

ADD_DEFINITIONS(-DUSEOMP)
ADD_DEFINITIONS(-DUSE_OPENMP)
ADD_DEFINITIONS(-DMULTI_THREADED)
ADD_DEFINITIONS(-DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS)
ADD_DEFINITIONS(-DBOOST_ERROR_CODE_HEADER_ONLY)
ADD_DEFINITIONS(-DBOOST_SYSTEM_NO_DEPRECATED)
ADD_DEFINITIONS(-DCXX11)
ADD_DEFINITIONS(-DPTR_STD)

SET_SOURCE_FILES_PROPERTIES( ../../src/protocols/cryst/spacegroup.cc PROPERTIES COMPILE_FLAGS -O0 )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -ffor-scope -std=c++11 -pipe -O0 -g")
SET(CMAKE_SHARED_LINKER_FLAGS "-lgomp")
SET(CMAKE_EXE_LINKER_FLAGS "-lgomp  -L/usr/local/lib")
include_directories("/usr/local/include")


INCLUDE(../build/build.post.cmake)

############ Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )

