CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)

#SET(COMPILER clang)
SET(COMPILER gcc)
SET(LINK_RULES shared)
#SET(LINK_RULES static)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)

SET(COMPILE_FLAGS "-std=c++98 -pipe -ffor-scope -w -O3 -ffast-math -funroll-loops")
#SET(COMPILE_FLAGS "-pipe -w")
#SET(COMPILE_FLAGS "-pipe -w -O2")
#SET(COMPILE_FLAGS "-pipe -w -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -w -ffast-math -funroll-loops -finline-functions")

#IF(${COMPILER} STREQUAL "clang")
#	SET(CMAKE_CXX_COMPILER "clang++")
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4")
#	#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#ELSE(${COMPILER} STREQUAL "clang")
#	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -finline-limit=20000 -s")
#ENDIF(${COMPILER} STREQUAL "clang")

INCLUDE(../build/build.post.cmake)

## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
