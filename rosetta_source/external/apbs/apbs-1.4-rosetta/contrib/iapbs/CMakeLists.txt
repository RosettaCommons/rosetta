################################################################################
# Set up basic project stuff                                                   #
################################################################################
cmake_minimum_required(VERSION 2.6)

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(CheckTypeSize)
include(ExternalProject)

set(iAPBS_VERSION "2.1.0")
#set(PACKAGE_STRING ${iAPBS_VERSION})

PROJECT(iAPBS)
#PROJECT(iAPBS C Fortran)
#ENABLE_LANGUAGE(Fortran)
################################################################################
# Set project paths                                                            #
################################################################################

message(STATUS "Setting iAPBS project paths")
set(iAPBS_ROOT ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${iAPBS_ROOT}/bin)
set(LIBRARY_OUTPUT_PATH ${iAPBS_ROOT}/lib)

set(LIBRARY_INSTALL_DIR lib)

################################################################################
# add src
################################################################################


add_subdirectory(src)

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
