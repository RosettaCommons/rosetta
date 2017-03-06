CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)
SET(COMPILER mpi)
SET(LINK_RULES shared)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DUSEMPI)

SET(COMPILE_FLAGS "-pipe -w -O3 -ffast-math -fno-finite-math-only -funroll-loops -finline-functions -finline-limit=20000 -s")

# from build_staticmpi ...
#SET(LINK_DYNAMIC_LIBS "-Wl,-Bdynamic -l c")
#SET(TRAILING_STATIC "-Wl,-Bstatic")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -w -ffast-math -fno-finite-math-only -funroll-loops -finline-functions")

IF(${COMPILER} STREQUAL "clang")
	SET(CMAKE_CXX_COMPILER "clang++")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ELSE(${COMPILER} STREQUAL "clang")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -finline-limit=20000 -s")
ENDIF(${COMPILER} STREQUAL "clang")


INCLUDE( ../build/build.post.cmake )
## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
