CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(MODE release)

#SET(COMPILER clang)
SET(COMPILER gcc)
SET(LINK_RULES shared)
#SET(LINK_RULES static)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DGL_GRAPHICS)

find_package( GLUT REQUIRED )
find_package( OpenGL REQUIRED )

INCLUDE_DIRECTORIES( /usr/X11R6/include )
LINK_DIRECTORIES( /usr/X11R6/lib )

# horrific hack -- and it doesn't even work. Are these flags not being used?
SET( CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS}  -framework GLUT -framework OpenGL -dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib")

SET(COMPILE_FLAGS "-pipe -w ")
#SET(COMPILE_FLAGS "-pipe -w -O2")
#SET(COMPILE_FLAGS "-pipe -w -O3 -ffast-math -funroll-loops -finline-functions -finline-limit=20000 -s")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -ffast-math -funroll-loops -finline-functions -ftemplate-depth=250 -stdlib=libstdc++ -mmacosx-version-min=10.6")

IF(${COMPILER} STREQUAL "clang")
	SET(CMAKE_CXX_COMPILER "clang++")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4")
ELSE(${COMPILER} STREQUAL "clang")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
ENDIF(${COMPILER} STREQUAL "clang")

INCLUDE(../build/build.post.cmake)

## Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )
