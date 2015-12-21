PROJECT(minirosetta)

# if you don't want the full compiler output, remove the following line
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# explicitly set this policy to suppress
# warning in cmake versions>=2.8.10.1
CMAKE_POLICY( SET CMP0011 NEW )

# this seems to fix some linking errors
SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cq <TARGET> <LINK_FLAGS> <OBJECTS>")
SET(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> q  <TARGET> <LINK_FLAGS> <OBJECTS>")

# include / link dirs
INCLUDE_DIRECTORIES(../..)
INCLUDE_DIRECTORIES(../../src)

# GENERAL EXTERNAL DIRECTORY SETTINGS

CMAKE_POLICY( SET CMP0015 NEW )

# We use the SYSTEM directive here to avoid printing warning in the external libraries

INCLUDE_DIRECTORIES(SYSTEM ../../external)
INCLUDE_DIRECTORIES(SYSTEM ../../external/include)
INCLUDE_DIRECTORIES(SYSTEM ../../external/dbio)

ADD_DEFINITIONS(-DBOOST_ERROR_CODE_HEADER_ONLY)
ADD_DEFINITIONS(-DBOOST_SYSTEM_NO_DEPRECATED)
INCLUDE_DIRECTORIES(SYSTEM ../../external/boost_1_55_0)

INCLUDE_DIRECTORIES(SYSTEM ../../external/cxxtest)

LINK_DIRECTORIES(../../external/boost_1_55_0)
LINK_DIRECTORIES(../../external/lib)

# Platform-specific includes
if(APPLE)
	INCLUDE_DIRECTORIES(../../src/platform/macos)
	ADD_DEFINITIONS(-DMAC)
ELSEIF(UNIX)
	INCLUDE_DIRECTORIES(../../src/platform/linux)
	ADD_DEFINITIONS(-DLINUX)
ELSEIF(WIN32)
	INCLUDE_DIRECTORIES(../../src/platform/windows)
	ADD_DEFINITIONS(-DWIN32)
ENDIF(APPLE)
