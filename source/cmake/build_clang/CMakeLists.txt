CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)

SET(COMPILER clang)
#SET(COMPILER gcc)
SET(LINK_RULES shared)
#SET(LINK_RULES static)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)

SET(COMPILE_FLAGS "-pipe -w")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -w -ffast-math -funroll-loops -finline-functions")

IF(${COMPILER} STREQUAL "clang")
    SET(CMAKE_C_COMPILER  "clang")
	SET(CMAKE_CXX_COMPILER "clang++")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++98 -D__extern_always_inline=inline")
ENDIF(${COMPILER} STREQUAL "clang")

INCLUDE(../build/build.post.cmake)

## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
