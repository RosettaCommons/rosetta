CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

##############################################################################
# This file only needs to specify which mode (i.e., release, debug), 
# extras (i.e., graphics), and compiler (i.e., gcc, clang) should be used. The
# appropriate settings for each mode, extras, and compiler should be 
# implemented in ../build/build.settings.cmake.
##############################################################################

############ Build Type, Mode, Link Rules, Compiler, Extras
SET(MODE debug)
SET(LINK_RULES shared)
SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)


############ Build Rules
INCLUDE(../build/compiler.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)
INCLUDE(../build/build.settings.cmake)
INCLUDE(../build/build.post.cmake)

############ Uncomment the following lines if you only wish to
############ build the rosetta_scripts app.
#INCLUDE( ../build/apps/rosetta_scripts.cmake )

############ THE FOLLOWING LINES ARE THE DEFAULT BEHAVIOUR.
############ ALL BINARIES, ALL PILOT APPS, AND ALL UNIT TESTS
############ WILL BE BUILT.
############ Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )

############ Unit Tests
INCLUDE(../build/build.test.cmake)
