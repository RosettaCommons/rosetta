cmake_minimum_required (VERSION 2.8)

################################################################################
# TODO:  Add licensing and authorship information
# TODO:  Test with FEtk build
# TODO:  Handle special mac dependencies
#        (configure.ac:1306)
################################################################################

################################################################################
# Future Enchancements
# --------------------
# [ ] Adding a profiling mode set by a flag
# [ ] Add functionality for creating rpm and deb packages
################################################################################


################################################################################
# Test platforms
# --------------
# [x] Ubuntu x86_64
# [x] Ubuntu i386
# [ ] Redhat 5
# [x] Redhat 6
# [x] Mac OSX
# [ ] Windows 7 x86_64
# [ ] Windows 7 i386
# [x] Windows 7 with Cygwin
# [ ] Windows 7 with Mingw
################################################################################



################################################################################
# Set up basic project stuff                                                   #
################################################################################

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(ExternalProject)

set(APBS_VERSION "1.4")
set(PACKAGE_STRING ${APBS_VERSION})
set(CMAKE_BUILD_TYPE "RELWITHDEBINFO")

project(apbs)

################################################################################
# Set project paths                                                            #
################################################################################

message(STATUS "Setting project paths")

set(APBS_ROOT ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${APBS_ROOT}/bin)
set(LIBRARY_OUTPUT_PATH ${APBS_ROOT}/lib)
set(TOOLS_PATH ${APBS_ROOT}/tools)
set(APBS_BINARY ${EXECUTABLE_OUTPUT_PATH}/apbs)

set(LIBRARY_INSTALL_PATH lib)
set(HEADER_INSTALL_PATH include/apbs)
set(EXECUTABLE_INSTALL_PATH bin)
set(SHARE_INSTALL_PATH share/apbs)

find_file(CONTRIB_PATH "contrib" PATHS "${APBS_ROOT}"
    DOC "The path to contributed modules for apbs"
    )



################################################################################
# Set up temporary files and directories                                       #
################################################################################

file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/temp)



################################################################################
# Set the lookup paths for external library dependencies                       #
################################################################################

message(STATUS "Setting lookup paths for headers and libraries")

set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}")
list(APPEND CMAKE_INCLUDE_PATH /usr/include)
list(APPEND CMAKE_INCLUDE_PATH /usr/local/include)

set(APBS_LIBS "")


################################################################################
# Enable ansi pedantic compiling                                               #
################################################################################

option(ENABLE_PEDANTIC "Enable the pedantic ansi compilation" OFF)

if(ENABLE_PEDANTIC)
    ADD_DEFINITIONS("-Wall -pedantic -ansi")
    message(STATUS "Pedantic compilation enabled")
endif()



################################################################################
# Determine Machine Epsilon values                                             #
################################################################################

OPTION(CHECK_EPSILON "" YES)
if(CHECK_EPSILON)
    message(STATUS "Computing machine epsilon values")
    try_run(
        FLOAT_EPSILON_COMPILED
        FLOAT_EPSILON_COMPUTED
        ${EXECUTABLE_OUTPUT_PATH}/temp
        ${APBS_ROOT}/src/.config/float_epsilon.c
        RUN_OUTPUT_VARIABLE FLOAT_EPSILON_OUTPUT
    )

    if(FLOAT_EPSILON_COMPUTED)
        message(STATUS "Floating point epsilon is ${FLOAT_EPSILON_OUTPUT}")
        set(FLOAT_EPSILON ${FLOAT_EPSILON_OUTPUT})
    else()
        message(FATAL_ERROR "Couldn't compute floating point machine epsilon")
    endif()

    try_run(
        DOUBLE_EPSILON_COMPILED
        DOUBLE_EPSILON_COMPUTED
        ${EXECUTABLE_OUTPUT_PATH}/temp
        ${APBS_ROOT}/src/.config/double_epsilon.c
        RUN_OUTPUT_VARIABLE DOUBLE_EPSILON_OUTPUT
    )

    if(DOUBLE_EPSILON_COMPUTED)
        message(STATUS "Double precision epsilon is ${DOUBLE_EPSILON_OUTPUT}")
        set(DOUBLE_EPSILON ${DOUBLE_EPSILON_OUTPUT})
    else()
        message(FATAL_ERROR "Couldn't compute double precision machine epsilon")
    endif()
else()
        set(FLOAT_EPSILON "1.19209290e-7")
        set(DOUBLE_EPSILON "2.2204460492503131e-16")
endif()

################################################################################
# Find some libraries                                                          #
################################################################################

if(NOT WIN32)
    find_library(MATH_LIBRARY "m")
    list(APPEND APBS_LIBS ${MATH_LIBRARY})
endif()



################################################################################
# Check for a few required functions                                           #
################################################################################

CHECK_FUNCTION_EXISTS(time HAVE_TIME_FUNC)

if(NOT HAVE_TIME_FUNC)
    message(FATAL_ERROR "Required time function not found")
endif()



CHECK_FUNCTION_EXISTS(rand HAVE_RAND_FUNC)

if(NOT HAVE_RAND_FUNC)
    message(FATAL_ERROR "Required rand function not found")
endif()



CHECK_FUNCTION_EXISTS(srand HAVE_SRAND_FUNC)

if(NOT HAVE_SRAND_FUNC)
    message(FATAL_ERROR "Required srand function not found")
endif()



################################################################################
# Look for the FEtk library                                                    #
################################################################################

find_file(
    FETK_PATH
    NAMES fetk
    PATHS ${SYS_LIBPATHS} /data/work/source
    DOC "Use an external fetk library"
    )

if(FETK_PATH)
    include_directories(${FETK_PATH}/include)
    message(STATUS "Found fetk: ${FETK_PATH}")
else()
    message(WARNING "Did not find fetk")
endif()



################################################################################
# Handle the MALOC library dependency                                          #
################################################################################

CHECK_INCLUDE_FILES(maloc/maloc.h HAVE_MALOC_H)

if(HAVE_MALOC_H)
    message(STATUS "External maloc headers found")
else()
    message(WARNING "External maloc headers not found.")
endif()

find_library(MALOC_LIBRARY  "maloc"
    PATHS "${FETK_PATH}/lib" "${CONTRIB_PATH}"
    DOC   "The fetk/maloc library"
)

if(MALOC_LIBRARY)
    message(STATUS "External maloc library was found: ${MALOC_LIBRARY}")
else()
    message(WARNING "External maloc library was not found.")
endif()

if(NOT (HAVE_MALOC_H AND MALOC_LIBRARY))
    message(STATUS "Using internal contributed maloc library")
    ExternalProject_add(
        "maloc"
        URL file://${CONTRIB_PATH}/maloc.tar.gz
        SOURCE_DIR ${CONTRIB_PATH}/maloc
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${APBS_ROOT}
-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
    )
    include_directories(${APBS_ROOT}/include)
    if(BUILD_SHARED_LIBS)
        set(MALOC_LIBRARY_BASENAME ${CMAKE_SHARED_LIBRARY_PREFIX}maloc${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
        set(MALOC_LIBRARY_BASENAME ${CMAKE_STATIC_LIBRARY_PREFIX}maloc${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()
    install(
        FILES ${LIBRARY_OUTPUT_PATH}/${MALOC_LIBRARY_BASENAME}
        DESTINATION ${LIBRARY_INSTALL_PATH}
    )
    set(MALOC_LIBRARY ${LIBRARY_OUTPUT_PATH}/${MALOC_LIBRARY_BASENAME})
endif()
message(STATUS "External maloc library was found: ${MALOC_LIBRARY}")
list(APPEND APBS_LIBS ${MALOC_LIBRARY})



################################################################################
# Handle the finite element solver dependencies                                #
################################################################################

option(ENABLE_FETK "Enable the finite element solver" OFF)

if(ENABLE_FETK)
    message(STATUS "Checking for fetk components")
    set(FETK_ENALBED 1)

    find_library(
        PUNC_LIBRARY
        NAMES "punk"
        PATHS ${FETK_PATH}/lib ${SYS_LIBPATHS}
        DOC   "The fetk/punc library"
        )

    if(PUNC_LIBRARY)
        list(APPEND APBS_LIBS ${PUNC_LIBRARY})
        set(HAVE_PUNC_H 1)
        message(STATUS "The punc library was found: ${PUNC_LIBRARY}")
    else()
        message(FATAL_ERROR "The punc library was not found")
    endif()



    find_library(
        MC_LIBRARY
        NAMES "mc"
        PATHS ${FETK_PATH}/lib ${SYS_LIBPATHS}
        DOC   "The fetk/mc library"
        )

    if(MC_LIBRARY)
        list(APPEND APBS_LIBS ${MC_LIBRARY})
        set(HAVE_MC_H 1)
        message(STATUS "The mc library was found: ${MC_LIBRARY}")
    else()
        message(FATAL_ERROR "The mc library was not found")
    endif()



    find_library(
        MCX_LIBRARY
        NAMES "mcx"
        PATHS ${FETK_PATH}/lib ${SYS_LIBPATHS}
        DOC   "The fetk/mcx library"
        )

    if(MCX_LIBRARY)
        list(APPEND APBS_LIBS ${MCX_LIBRARY})
        set(HAVE_MCX_H 1)
        message(STATUS "The mcx library was found: ${MCX_LIBRARY}")
    else()
        message(FATAL_ERROR "The mcx library was not found")
    endif()



    find_library(
        GAMER_LIBRARY
        NAMES "gamer"
        PATHS ${FETK_PATH}/lib ${SYS_LIBPATHS}
        DOC   "The fetk/gamer library"
        )

    if(GAMER_LIBRARY)
        list(APPEND APBS_LIBS ${GAMER_LIBRARY})
        set(HAVE_BIOM_H 1)
        message(STATUS "The gamer library was found: ${GAMER_LIBRARY}")
    else()
        message(FATAL_ERROR "The gamer library was not found")
    endif()



    find_library(
        SUPERLU_LIBRARY
        NAMES "superlu"
        PATHS ${SYS_LIBPATHS}
        DOC "The superlu library"
    )

    if(SUPERLU_LIBRARY)
        list(APPEND APBS_LIBS ${SUPERLU_LIBRARY})
    else()
        message(FATAL_ERROR "The superlu library was not found")
    endif()



    find_library(
        UMFPACK_LIBRARY
        NAMES "umfpack"
        PATHS ${SYS_LIBPATHS}
        DOC "The umfpack library"
    )

    if(UMFPACK_LIBRARY)
        list(APPEND APBS_LIBS ${UMFPACK_LIBRARY})
    else()
        message(FATAL_ERROR "The umfpack library was not found")
    endif()



    find_library(
        BLAS_LIBRARY
        NAMES "blas"
        PATHS ${SYS_LIBPATHS}
        DOC "The blas library"
    )

    if(BLAS_LIBRARY)
        list(APPEND APBS_LIBS ${BLAS_LIBRARY})
    else()
        message(FATAL_ERROR "The blas library was not found")
    endif()

endif()


################################################################################
# Handle conditional fast mode                                                 #
################################################################################

option(ENABLE_FAST "Enable fast mode" OFF)

if(ENABLE_FAST)
    set(APBS_FAST 1)
    message(STATUS "Fast mode enabled")
endif()



################################################################################
# Handle conditional TINKER support                                            #
################################################################################

option(ENABLE_TINKER "Enable TINKER support" OFF)

if(ENABLE_TINKER)
    set(WITH_TINKER 1)
    message(STATUS "Tinker enabled")
endif()



################################################################################
# Handle conditional availability of macro embedding                           #
################################################################################

try_compile(
    HAVE_EMBED
    ${APBS_ROOT}/build
    ${APBS_ROOT}/src/.config/embed_test.c
)

# TODO: Determine if the EMBED macro is even used



################################################################################
# Handle conditional debug building                                            #
################################################################################

option(ENABLE_DEBUG "Enable debugging mode" OFF)

if(ENABLE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    set(DEBUG 1)
    message(STATUS "Debugging compilation enabled")
endif()



################################################################################
# Enable inline functions conditionally dependant on debug mode                #
################################################################################


option(ENABLE_INLINE "Enable inline functions" ON)

if(ENABLE_INLINE)
    set(APBS_NOINLINE 1)
    message(STATUS "Inline functions enabled")
endif()



################################################################################
# Handle conditional building with quiet mode                                  #
################################################################################

option(ENABLE_QUIET "Enable quiet mode" OFF)

if(ENABLE_QUIET)
    set(VAPBSQUIET 1)
    message(STATUS "Quiet mode enabled")
endif()



################################################################################
# Handle conditional building with verbose debugging information printouts     #
################################################################################

option(ENABLE_VERBOSE_DEBUG "Enable verbose debugging mode" OFF)

if(ENABLE_VERBOSE_DEBUG)
    set(VERBOSE_DEBUG 1)
    message(STATUS "Verbose debugging mode enabled")
endif()



################################################################################
# Configure Python Wrappers                                                    #
################################################################################

option(ENABLE_PYTHON "Enable python wrappers" OFF)
option(MAX_MEMORY "Set the maximum memory (in MB) to be used for a job" OFF)

if(ENABLE_PYTHON)

    message(STATUS "Looking for python components")

    if(NOT CMAKE_COMPILER_IS_GNUCC)
        message(FATAL_ERROR "Python wrappers require GNU compilers")
    endif()


    find_package(PythonInterp)
    if(NOT PYTHON_EXECUTABLE)
        message(FATAL_ERROR "Couldn't find python interpreter")
    endif()

    find_package(PythonLibs)
    if(NOT PYTHON_LIBRARY)
        message(FATAL_ERROR "Couldn't find python libraries")
    endif()

    # TODO: Handle special flags for Python linking on different platforms
    # TODO: Handle special flags for Python linking with intel compiler
    # ${PYTHON_LFLAGS} perhaps?

    set(MODULE_PROBE_PY "")
    list(APPEND MODULE_PROBE_PY "from sys import stdout")
    list(APPEND MODULE_PROBE_PY "from distutils.sysconfig import get_python_lib")
    list(APPEND MODULE_PROBE_PY "stdout.write( get_python_lib() )")
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "${MODULE_PROBE_PY}"
        ERROR_QUIET
        OUTPUT_VARIABLE PYTHON_MODULES
    )

    if(NOT PYTHON_MODULES)
        message(FATAL_ERROR "Coudln't determine location of python modules")
    endif()

    if(NOT MAX_MEMORY)
        set(MAX_MEMORY "-1")
    endif()


    # TODO:  these can't e right...
    set(SERVICEURL "http://kryptonite.nbcr.net/opal2/services/apbs_1.3")
    set(PARALLELSERVICEURL "http://oolite.calit2.optiputer.net/opal2/services/apbs-parallel-1.3")



    configure_file(
        ${APBS_ROOT}/src/.config/ApbsClient.py.in
        ${EXECUTABLE_OUTPUT_PATH}/ApbsClient.py
    )

endif()



################################################################################
# Handle conditional building with OpenMP                                      #
################################################################################

option(ENABLE_OPENMP "Enable OpenMP parallelism" ON)

if(ENABLE_OPENMP)
    if(NOT ENABLE_DEBUG)
        message(STATUS "Checking for OpenMP")
        find_package(OpenMP)
        if(OPENMP_FOUND)
            message(STATUS "OpenMP support enabled")
            add_definitions("${OpenMP_C_FLAGS}")
            list(APPEND APBS_LIBS ${OpenMP_C_FLAGS})
        else()
            message(WARNING "OpenMP was not found.  OpenMP support disabled")
        endif()
    else()
        message(WARNING "OpenMP may not be enabled in debugging mode")
    endif()
endif()



################################################################################
# Handle conditional building with MPI Support                                 #
################################################################################

option(ENABLE_MPI "Enable MPI parallelism" OFF)

if(ENABLE_MPI)
    if(NOT ENABLE_DEBUG)
        message(STATUS "Checking for MPI")
        find_package(MPI)
        if(MPI_FOUND)
            message(STATUS "MPI support enabled")
            include_directories(${MPI_INCLUDE_PATH})
            list(APPEND APBS_LIBS ${MPI_LIBRARIES})
        else()
            message(WARNING "MPI was not found; disabling")
        endif()
    else()
        message(WARNING "MPI may not be enabled in debugging mode")
    endif()
endif()



################################################################################
# Handle library checks for embedded unix environments in windows              #
################################################################################

if(MINGW)
    message(STATUS "Checking for wsock32 in MinGW environment")
    find_library(
        MINGW_WSOCK32
        NAMES wsock32
        PATHS ${SYS_LIBPATHS}
        DOC   "The wsock32 library"
        )

    if(MINGW_WSOCK32)
        message(STATUS "The wsock32 library was found: ${MINGW_HAS_WSOCK32}")
    else()
        message(FATAL_ERROR "The wsock32 library was not fond")
    endif()
endif()



if(CYGWIN)
    message(STATUS "Checking for wsock32 in Cygwin environment")
    find_library(
        CYGWIN_WSOCK32
        NAMES wsock32
        PATHS ${SYS_LIBPATHS}
        DOC   "The wsock32 library"
        )

    if(CYGWIN_WSOCK32)
        message(STATUS "The wsock32 library was found: ${CYGWIN_WSOCK32}")
        list(APPEND APBS_LIBS ${CYGWIN_WSOCK32})
    else()
        message(FATAL_ERROR "The wsock32 library was not fond")
    endif()

    set(HAVE_CYGWIN 1)
endif()

if(NOT CYGWIN AND NOT MINGW AND WIN32)
    list(APPEND APBS_LIBS wsock32)
    ADD_DEFINITIONS("/D _CRT_SECURE_NO_WARNINGS")
endif()



################################################################################
# Build APBS sources                                                           #
################################################################################

include_directories(${APBS_ROOT}/src)
add_subdirectory(src)



################################################################################
# Build APBS documentation                                                     #
################################################################################

option(BUILD_DOC "Build/Rebuild documentation using doxygen" OFF)

if(BUILD_DOC)
    message(STATUS "Building documentation enabled")
    add_subdirectory(doc)
endif()



################################################################################
# Handle conditional building with verbose debugging information printouts     #
################################################################################

option(BUILD_TOOLS "Build supplemental tools" ON)

if(BUILD_TOOLS)
    message(STATUS "Supplemental tools enabled")
    add_subdirectory(tools)
endif()



################################################################################
# Set up additional directories to install                                     #
################################################################################

install(
    DIRECTORY ${APBS_ROOT}/doc
    DESTINATION ${SHARE_INSTALL_PATH}
)

install(
    DIRECTORY ${APBS_ROOT}/tests
    DESTINATION ${SHARE_INSTALL_PATH}
    FILES_MATCHING
        PATTERN "*.py"
        PATTERN "README"
)

if(NOT HAVE_MALOC_H)
    message(STATUS "Installing maloc headers")
    install(
        DIRECTORY ${APBS_ROOT}/include/maloc
        DESTINATION ${HEADER_INSTALL_PATH}/..
    )
endif()

################################################################################
# Optionally build iAPBS interface                                             #
################################################################################

option(ENABLE_iAPBS "Enable iAPBS" OFF)

if(ENABLE_iAPBS)
    message(STATUS "Building of iAPBS interface enabled")
    add_subdirectory(contrib/iapbs/src)
    list(APPEND APBS_LIBS ${iAPBS_LIBRARY})
endif()



################################################################################
# Clean up temporary files and directories                                     #
################################################################################

file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}/temp)
