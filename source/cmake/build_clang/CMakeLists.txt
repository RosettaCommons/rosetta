CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)

SET(COMPILER clang)
#SET(COMPILER gcc)
SET(LINK_RULES shared)
#SET(LINK_RULES static)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)

SET(COMPILE_FLAGS "-pipe -w -O3 -std=c++98 -Qunused-arguments -D__extern_always_inline=inline")

INCLUDE(../build/build.post.cmake)

## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
