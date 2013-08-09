CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE debug)
SET(COMPILER gcc)
#SET(COMPILER clang)
SET(LINK_RULES shared)
#SET(LINK_RULES static)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-D_DEBUG)
SET(COMPILE_FLAGS "-c -std=c++98 -pipe -ffor-scope -Wall -Wextra -pedantic -Wno-long-long -O0 -g -ggdb -ffloat-store")

INCLUDE(../build/build.post.cmake)

## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
