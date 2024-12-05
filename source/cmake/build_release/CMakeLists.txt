CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

##############################################################################
# This file only needs to specify which mode (i.e., release, debug), 
# extras (i.e., graphics), and compiler (i.e., gcc, clang) should be used. The
# appropriate settings for each mode, extras, and compiler should be 
# implemented in ../build/build.settings.cmake.
##############################################################################

############ Build Type, Mode, Link Rules, Compiler, Extras
SET(MODE release)
SET(LINK_RULES shared)
SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)

############ Build Rules
INCLUDE(../build/compiler.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/${LINK_RULES}.build.cmake)
INCLUDE(../build/build.settings.cmake)
INCLUDE(../build/build.post.cmake)

########### Uncomment the following if you only want to build RosettaScripts:
#INCLUDE( ../build/apps/rosetta_scripts.cmake )

############ Top level binaries
INCLUDE( ../build/apps.all.cmake )
INCLUDE( ../build/pilot_apps.all.cmake )

############ Unit Tests
INCLUDE(../build/build.test.cmake)

############ Installation rules
INCLUDE( ../build/install.cmake )
