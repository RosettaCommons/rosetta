# The settings below are combined together to create the final set of settings
# for CMake + Ninja builds. These settings attempt to match those found in
# tools/build/basic.settings, which are used for building with scons.
#
# Currently, this file only includes settings for CMake + Ninja release, debug,
# and graphics builds. This file is included in the CMakeLists.txt files found
# in the build_release, build_debug, and build_graphics directories. The
# CMakeLists.txt files only need to specify which mode (i.e., release, debug),
# extras (i.e., graphics), and compiler (i.e., gcc, clang) should be used. The
# appropriate settings for each mode, extras, and compiler should be
# implemented here.
#
# Ultimately, it may be more efficient to automate the generation of this
# build.settings.cmake file so that it does not have to be updated manually
# every time basic.settings is modified.
#
# @file source/cmake/build/build.settings.cmake
# @author Caleb W. Geniesse (calebgeniesse@stanford.edu)


###########################################################################
# GCC #####################################################################
###########################################################################

if( ${COMPILER} STREQUAL "gcc" )
	ADD_DEFINITIONS(-DPTR_STD)
	ADD_DEFINITIONS(-DCXX11)
	set( cc
			#-std=c99
	)
	set( cxx
			-std=c++11
	)
	set( compile
			-pipe
			-ffor-scope
			-ftemplate-depth-256
			-fPIC
			-DBOOST_ERROR_CODE_HEADER_ONLY
			-DBOOST_SYSTEM_NO_DEPRECATED
			-I /usr/include
			-I /usr/local/include
			-I src
			-I external/include
			-I src/platform/linux
	)
	set( warn
			-Wall
			-Wextra
			-pedantic
			-Werror  # REMOVE FOR RELEASE
			-Wno-long-long
			-Wno-strict-aliasing
	)

endif()

# modes ###################################################################

# "gcc, debug"
if( ${COMPILER} STREQUAL "gcc" AND ${MODE} STREQUAL "debug" )
	list( APPEND compile
			-O0
	)
	list( APPEND mode
			-g
			-ggdb
			-ffloat-store
			#-fstack-check
	)
	list( APPEND defines
			-D_DEBUG
	)
endif()

# "gcc, release"
if( ${COMPILER} STREQUAL "gcc" AND ${MODE} STREQUAL "release" )
	list( APPEND compile
			-O3
			-ffast-math
			-fno-finite-math-only
			-funroll-loops
			-finline-functions
			-finline-limit=20000
			-s
	)
	list( APPEND warn
			-Wno-unused-variable
			-Wno-unused-parameter
	)
	list( APPEND defines
			-DNDEBUG
	)
endif()

# "gcc", "release_bluegene"
# Used exclusively for compilation on the Argonne "Mira" Blue Gene supercomputer.
# Added by Vikram K. Mulligan, Baker lab (vmullig@uw.edu) on 19 April 2016.
if( ${COMPILER} STREQUAL "gcc" AND ${MODE} STREQUAL "release_bluegene" )
	list( REMOVE_ITEM compile
			-I /usr/include
			-I /usr/local/include
			-I src
			-I external/include
			-I src/platform/linux

	)
	list( APPEND compile
			-I ../../external/boost_1_55_0
			-I ../../external/include
			-I ../../src
			-I ../../src/platform/linux
			#-I /usr/include
			#-I /usr/local/include
			-O3
			-ffast-math
			-fno-finite-math-only
			-funroll-loops
			-finline-functions
			-finline-limit=20000
			-s
	)
	list ( APPEND warn
			-Wno-deprecated
			-Wno-unused-variable
			-Wno-unused-parameter
			-Wno-type-limits
	)
	list( APPEND defines
			-DNDEBUG
			#-DDISABLE_SQLITE
	)
endif()

###########################################################################
# Windows Visual Studio - CL ##############################################
###########################################################################

if( ${COMPILER} STREQUAL "cl" )
	#list( REMOVE_ITEM warn
	#		-Wextra
	#)

endif()

###########################################################################
# Blue Gene XLC compiler ##################################################
###########################################################################

# "xlc", "release_bluegene"
# Used exclusively for compilation on the Argonne "Mira" Blue Gene supercomputer.
# Added by Vikram K. Mulligan, Baker lab (vmullig@uw.edu) on 20 April 2016.
if( ${COMPILER} STREQUAL "xlc" AND ${MODE} STREQUAL "release_bluegene" )
        ADD_DEFINITIONS(-DPTR_BOOST)
        set( cc
                        #-std=c99
        )
        set( cxx
                        -std=c++11
        )
        set( compile
			-O3
			-s
			#-DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
			#-DBOOST_ERROR_CODE_HEADER_ONLY
			-DBOOST_SYSTEM_NO_DEPRECATED
			-ftemplate-depth-256
			-fPIC
			-DBOOST_ERROR_CODE_HEADER_ONLY
			-DBOOST_SYSTEM_NO_DEPRECATED
			-I /usr/include
			-I /usr/local/include
			-I src
			-I external/include
			-I src/platform/linux
        )
        list( APPEND defines
                        -DNDEBUG
        )
endif()

###########################################################################
# Clang ###################################################################
###########################################################################

set(WITH_LIBSTDC++ ON CACHE BOOL "Build libraries using libstdc++ when using clang compiler. ")

if( ${COMPILER} STREQUAL "clang" )
	ADD_DEFINITIONS(-DPTR_STD)
	ADD_DEFINITIONS(-DCXX11)
	set( cc
			#-std=c99
	)
	set( cxx
			-std=c++11
	)
	set( compile
			-pipe
			-Qunused-arguments
			-DUNUSUAL_ALLOCATOR_DECLARATION
			-ftemplate-depth-256
			-stdlib=libstdc++
			-fPIC
			-DBOOST_ERROR_CODE_HEADER_ONLY
			-DBOOST_SYSTEM_NO_DEPRECATED
			-I /usr/include
			-I /usr/local/include
			-I src
			-I external/include
			-I src/platform/linux
			-I src/platform/macos/64/clang/6.1
			-I src/platform/macos/64/clang
			-I src/platform/macos/64
	)
	set( warn
			-W
			-Wall
			-Wextra
			-pedantic
			#-Weverything
			-Werror # REMOVE FOR RELEASE
			-Wno-long-long
			-Wno-strict-aliasing
			#-Wno-documentation
			#-Wno-padded
			#-Wno-weak-vtables
	)


	set( shlink
		-stdlib=libstdc++
	)

	set( link
		-stdlib=libstdc++
	)

	if( WITH_LIBSTDC++ )
		list( APPEND compile
			-stdlib=libstdc++
		)
		list (APPEND shlink
			-stdlib=libstdc++
		)
		list (APPEND link
			-stdlib=libstdc++
		)
	endif()
endif()

if(APPLE AND ${COMPILER} STREQUAL "clang")
    list( REMOVE_ITEM compile
	-stdlib=libstdc++
    )
    list( REMOVE_ITEM shlink
	-stdlib=libstdc++
    )
    list( REMOVE_ITEM link
	-stdlib=libstdc++
    )

    set( cxx
          -std=c++11
	  -stdlib=libc++
    )
    set( shlink
	  -stdlib=libc++
    )
    set( link
	  -stdlib=libc++
    )

endif()


# modes ###################################################################

# "clang, debug"
if( ${COMPILER} STREQUAL "clang" AND ${MODE} STREQUAL "debug" )
	list( APPEND compile
			-O0
	)
	list( APPEND mode
			-g
	)
endif()

# "clang, release"
if( ${COMPILER} STREQUAL "clang" AND ${MODE} STREQUAL "release" )
	list( APPEND compile
			-O3
	)
	list( APPEND warn
			-Wno-unused-variable
			-Wno-unused-parameter
	)
	list( APPEND defines
			-DNDEBUG
	)

endif()

# "clang", "release_bluegene"
# Used exclusively for compilation on the Argonne "Mira" Blue Gene supercomputer.
# Added by Vikram K. Mulligan, Baker lab (vmullig@uw.edu) on 11 Sept 2016.
if( ${COMPILER} STREQUAL "clang" AND ${MODE} STREQUAL "release_bluegene" )
	list( REMOVE_ITEM compile
			-I /usr/include
			-I /usr/local/include
			-I src
			-I external/include
			-I src/platform/linux
			-I src/platform/macos/64/clang/6.1
			-I src/platform/macos/64/clang
			-I src/platform/macos/64
			-stdlib=libstdc++
	)
	list( APPEND compile
			-I ../../external/boost_1_55_0
			-I ../../external/include
			-I ../../src
			#-I /usr/include
			#-I /usr/local/include
			-I ../../src/platform/linux
			-O3
			-stdlib=libc++
	)
	list( REMOVE_ITEM cxx
			-stdlib=libstdc++
	)
	set( cxx
			-std=c++11
			-stdlib=libc++
	)
	list( REMOVE_ITEM link
			-stdlib=libstdc++
	)
	list( APPEND link
			-lpthread
			-stdlib=libc++
	)
	list( REMOVE_ITEM shlink
			-stdlib=libstdc++
	)
	list( APPEND shlink
			-lpthread
			-stdlib=libc++
	)
	list ( APPEND warn
			-Wno-tautological-constant-out-of-range-compare
			-Wno-undefined-var-template
			-Wno-unused-variable
			-Wno-unused-parameter
			-Wno-type-safety
			-Wno-inconsistent-missing-override
	)
	list( APPEND defines
			-DNDEBUG
			-DBLUEGENECLANG
			#-DDISABLE_SQLITE
	)
endif()


# Pyrosetta build options
set(WITH_PYROSETTA OFF CACHE BOOL "Build rosetta libraries with pyrosetta-specific build configuration.")
if( WITH_PYROSETTA )
  list( APPEND defines -DPYROSETTA )
endif()

###########################################################################
# Extras ###################################################################
###########################################################################

if( EXTRAS )

	# "macos, graphics"
	if( APPLE AND ${EXTRAS} STREQUAL "graphics" )
		#find_package( GLUT REQUIRED )
		#find_package( OpenGL REQUIRED )
		if (EXISTS /usr/X11R6/include )
			include_directories( /usr/X11R6/include )
			link_directories( /usr/X11R6/lib )
		else()
			include_directories( /usr/X11/include )
			link_directories( /usr/X11/lib )
		endif()
		list( APPEND defines
				-DGL_GRAPHICS
				-DMAC
		)
		set( link
				-stdlib=libstdc++
				-stdlib=libc++
				-framework GLUT
				-framework OpenGL
				-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib
		)
		set( shlink
				-stdlib=libstdc++
				-stdlib=libc++
				-framework GLUT
				-framework OpenGL
				-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib
		)
		list( APPEND compile
				-pipe
				-ffast-math
				-fno-finite-math-only
				-funroll-loops
				-ftemplate-depth=250
				-mmacosx-version-min=10.10
		)
		set( warn
#				-Werror=sign-compare
				-Werror=reorder
				-Werror=address
				-Werror=char-subscripts
				-Werror=comment
				-Werror=nonnull
				-Werror=type-limits
				-Werror=parentheses
				-Werror=ignored-qualifiers
				-Werror=enum-compare
		)
	endif()

	# "linux, graphics"
	if( UNIX AND NOT APPLE AND ${EXTRAS} STREQUAL "graphics" )
		list( APPEND defines
				-DGL_GRAPHICS
		)
		list( APPEND link
				-lpthread
				-lGL
				-lGLU
				-lglut
		)
		list( APPEND shlink
				-lpthread
				-lGL
				-lGLU
				-lglut
		)
		set( warn
				-Werror=sign-compare
				-Werror=reorder
				-Werror=address
				-Werror=char-subscripts
				-Werror=comment
				-Werror=nonnull
				-Werror=type-limits
				-Werror=parentheses
				-Werror=ignored-qualifiers
				-Werror=enum-compare
		)
	endif()

endif()


###########################################################################
###########################################################################
###########################################################################

add_definitions( ${defines} )

foreach( flag ${shlink} )
	set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${flag}" )
endforeach()

foreach( flag ${cc} ${cxx} ${compile} ${warn} ${mode} )
	set( COMPILE_FLAGS "${COMPILE_FLAGS} ${flag}" )
endforeach()

foreach( flag ${COMPILE_FLAGS} )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" )
endforeach()
