# Require a fairly recent version of cmake.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Tell cmake where additional modules can be found.
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../build/finders)

# Specify the default compile flags.  Different options may add new flags to 
# this list.

SET(COMPILE_FLAGS "${COMPILE_FLAGS} -pipe -w -O3 -s")
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -ffast-math -fno-finite-math-only -funroll-loops")
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -finline-limit=20000")
SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Werror")

# Specify which compiler to use.  This option must be set before PROJECT() is 
# called, which happens on the first INCLUDE().  Note that it's considered bad 
# form to set the compiler like this, because it makes the build dependent on 
# 'clang++' being defined the same on everyone's computer.  Better is to use a 
# command-line argument when calling cmake.

SET(CMAKE_CXX_COMPILER "clang++")

# Set defaults for the variables that will be used to name the executables.  
# Note that these variables only affect the name of the generated executables, 
# even though they're named like they would also affect the build. 

SET(MODE "default_mode")
SET(COMPILER "default_compiler")
SET(LINK_RULES "shared")

# Setup custom build options.  Feel free to add new build options, but follow 
# the pattern below and create a new script for each option to keep this file 
# clean.  Some useful commands include:
#
# 1. SET(COMPILE_FLAGS) to new compilation flags.
# 2. INCLUDE_DIRECTORIES() to add external headers.
# 3. SET(LINK_EXTERNAL_LIBS) to add external libraries.
# 4. ADD_DEFINITIONS() to add preprocessor #defines.

INCLUDE(modules/compiler.cmake)
INCLUDE(modules/build_type.cmake)

# Define how both core rosetta libraries and custom build modules should be 
# built and linked.

INCLUDE(build.pre.cmake)
INCLUDE(${LINK_RULES}.build.cmake)
INCLUDE(modules/mpi.cmake)
INCLUDE(modules/mysql.cmake)
INCLUDE(modules/hdf5.cmake)
INCLUDE(build.post.cmake)

# Build the top-level binaries.

INCLUDE(apps.all.cmake)
INCLUDE(pilot_apps.all.cmake)

