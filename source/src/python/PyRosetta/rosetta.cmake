cmake_minimum_required(VERSION 3.0)

# enabling IN_LIST operator https://cmake.org/cmake/help/v3.7/policy/CMP0057.html
# cmake_policy(SET CMP0057 NEW)

# string below should form comment unique for given build config (compiler, gcc-install-prefix, python-include-dir, python-lib) options (we use this for dependency tracking
#%__PyRosetta_build_config__%#

# Export compilation database for use in IDEs
SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)

project(rosetta)

# Add a CMake parameter for choosing a desired Python version
set(PYROSETTA_PYTHON_VERSION "" CACHE STRING "Python version to use for compiling the PyRosetta library")
set(PYROSETTA_STRIP_MODULE TRUE CACHE BOOL "Strip compiled modules in release mode to minimize file size.")

include(CheckCXXCompilerFlag)

# Set a default build configuration if none is specified. 'MinSizeRel' produces the smallest binaries
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'MinSizeRel' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
# if( CMAKE_BUILD_TYPE MATCHES "Debug" OR CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo") # if( CMAKE_BUILD_TYPE IN_LIST "Debug" "RelWithDebInfo")
#   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcmodel=large")
#   add_compile_options(-mcmodel=large)
# endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

# Try to autodetect Python (can be overridden manually if needed)
if (NOT ${PYROSETTA_PYTHON_VERSION} STREQUAL "")
  find_package(PythonLibs ${PYROSETTA_PYTHON_VERSION} EXACT)
  if (NOT PythonLibs_FOUND)
     MESSAGE( FATAL_ERROR "Could not find requested Python version EXACT: " ${PYROSETTA_PYTHON_VERSION} " terminating...")
     find_package(PythonLibs ${PYROSETTA_PYTHON_VERSION} REQUIRED)
  endif()
else()
    set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6 3.7 2.7)
    find_package(PythonLibs REQUIRED)
endif()

# Locate zlib
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
     MESSAGE( FATAL_ERROR "Could not find zlib. terminating...")
endif()

# Uncomment the following line if you will also require a matching Python interpreter
# find_package(PythonInterp ${PYTHONLIBS_VERSION_STRING} EXACT REQUIRED)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_CPP11_FLAG)

  if (HAS_CPP14_FLAG  AND  APPLE)  # disabling C++14 support on Linux for now due to problem with GCC-4.8 libstdc++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  elseif (HAS_CPP11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else()
    message(FATAL_ERROR "Unsupported compiler -- at least C++11 support is needed!")
  endif()

  # Enable link time optimization and set the default symbol
  # visibility to hidden (very important to obtain small binaries)
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    # Default symbol visibility
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fvisibility=hidden -fvisibility-inlines-hidden")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden")


    # Check for Link Time Optimization support, for now disable this because it lead to multi-hour linking time on Mac
    # CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
    # if (HAS_LTO_FLAG)
    #   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    # endif()
  endif()
endif()

# if ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "9.0") )
# else()
#     add_compile_options(-Werror=non-virtual-dtor -Werror=delete-abstract-non-virtual-dtor) # -pedantic-errors -pedantic -Werror
# endif()

add_compile_options(-Werror=non-virtual-dtor) # -pedantic-errors -pedantic -Werror -Werror=delete-abstract-non-virtual-dtor


# work-around for: pybind11.h:147:36: error: invalid conversion from 'std::enable_if<true, void*>::type {aka void*}' to 'const pybind11::detail::void_type*' [-fpermissive]
# only works for GCC
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

# forcing  -fPIC so we can make a shared object
set_property(GLOBAL PROPERTY POSITION_INDEPENDENT_CODE ON)

# Include path for Python header files
include_directories(${PYTHON_INCLUDE_DIR})

# Include path for pybind11 header files -- this may need to be changed depending on your setup
#include_directories(${PROJECT_SOURCE_DIR}/pybind11/include)

#%__Rosetta_cmake_instructions__%#

# Code below moved to build.py
# Copy source and database files into build from source directory.
# file(COPY
#   src/python/PyRosetta/src/pyrosetta
#   src/python/PyRosetta/src/rosetta
#   src/python/PyRosetta/src/test
#   src/python/PyRosetta/src/demo
#   src/python/PyRosetta/src/self-test.py
#   src/python/PyRosetta/src/setup.py
#   src/python/PyRosetta/src/setup.cfg
#   src/python/PyRosetta/src/ez_setup.py
#   DESTINATION .
# )
# file(GLOB PYROSETTA_DATABASE_FILES database/*)
# file(COPY ${PYROSETTA_DATABASE_FILES} DESTINATION pyrosetta/database)

if(PYROSETTA_EXTERNAL_LINK)
    # Copy compiled external modules into pyrosetta module
    file(GLOB PYROSETTA_LIB_FILES lib/*.so)
    file(COPY ${PYROSETTA_LIB_FILES} DESTINATION pyrosetta/lib)
endif()

# Create the binding library
add_library(rosetta SHARED
#%__PyRosetta_sources__%#
)

target_compile_options(rosetta PRIVATE
  #%__PyRosetta_compile_options__%#
)

target_link_libraries(rosetta
  #%__PyRosetta_target_link_libraries__%#
)


target_link_libraries(rosetta #%__Rosetta_libraries__%#  ${ZLIB_LIBRARIES} )

# Don't add a 'lib' prefix to the shared library
set_target_properties(rosetta PROPERTIES PREFIX "")

if(PYROSETTA_EXTERNAL_LINK)
    #Need to set rpath property directly w/ link flag, as opposed to cmake rpath handling
    #as $ORIGIN is mangled when passed in via .rsp file when linking large libraries
    set_property(TARGET rosetta APPEND PROPERTY LINK_FLAGS "-Wl,-rpath,\\$ORIGIN/lib")
    #set_target_properties(rosetta  PROPERTIES INSTALL_RPATH "$ORIGIN")
    #set_target_properties(rosetta  PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
endif()

# Place compiled module into pyrosetta module directory
set_target_properties( rosetta
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyrosetta"
)

if (WIN32)
  if (MSVC)
    # /bigobj is needed for bigger binding projects due to the limit to 64k
    # addressable sections. /MP enables multithreaded builds (relevant when
    # there are many files).
    set_target_properties(rosetta PROPERTIES COMPILE_FLAGS "/MP /bigobj ")

    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      # Enforce size-based optimization and link time code generation on MSVC
      # (~30% smaller binaries in experiments).
      set_target_properties(rosetta APPEND_STRING PROPERTY COMPILE_FLAGS "/Os /GL ")
      set_target_properties(rosetta APPEND_STRING PROPERTY LINK_FLAGS "/LTCG ")
    endif()
  endif()

  # .PYD file extension on Windows
  set_target_properties(rosetta PROPERTIES SUFFIX ".pyd")

  # Link against the Python shared library
  target_link_libraries(rosetta ${PYTHON_LIBRARY})
elseif (UNIX)
  # It's quite common to have multiple copies of the same Python version
  # installed on one's system. E.g.: one copy from the OS and another copy
  # that's statically linked into an application like Blender or Maya.
  # If we link our plugin library against the OS Python here and import it
  # into Blender or Maya later on, this will cause segfaults when multiple
  # conflicting Python instances are active at the same time (even when they
  # are of the same version).

  # Windows is not affected by this issue since it handles DLL imports
  # differently. The solution for Linux and Mac OS is simple: we just don't
  # link against the Python library. The resulting shared library will have
  # missing symbols, but that's perfectly fine -- they will be resolved at
  # import time.

  # .SO file extension on Linux/Mac OS
  set_target_properties(rosetta PROPERTIES SUFFIX ".so")

  # Strip unnecessary sections of the binary on Linux/Mac OS
  if(APPLE)
    set_target_properties(rosetta PROPERTIES MACOSX_RPATH ".")
    set_target_properties(rosetta PROPERTIES LINK_FLAGS "-undefined dynamic_lookup ")
  endif()

  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG AND PYROSETTA_STRIP_MODULE)
    if(APPLE)
      add_custom_command(TARGET rosetta POST_BUILD COMMAND strip -u -r ${PROJECT_BINARY_DIR}/pyrosetta/rosetta.so)
    else()
      add_custom_command(TARGET rosetta POST_BUILD COMMAND strip ${PROJECT_BINARY_DIR}/pyrosetta/rosetta.so)
    endif()
  endif()

endif()
