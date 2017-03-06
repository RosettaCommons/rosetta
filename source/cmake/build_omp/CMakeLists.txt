CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)
SET(COMPILER gcc)
SET(LINK_RULES shared)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DUSEOMP)
ADD_DEFINITIONS(-DUSE_OPENMP)
ADD_DEFINITIONS(-DMULTI_THREADED)

SET(COMPILE_FLAGS "-pipe -w -O3 -ffast-math -fno-finite-math-only -funroll-loops -finline-functions -finline-limit=20000 -fopenmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -finline-limit=20000")
SET(CMAKE_SHARED_LINKER_FLAGS "-lgomp")
SET(CMAKE_EXE_LINKER_FLAGS "-lgomp")


INCLUDE( ../build/build.post.cmake )
## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
