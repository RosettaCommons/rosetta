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
	set( cc
			#-std=c99
			-isystem external/boost_1_55_0/
			-isystem external/include/
			-isystem external/dbio/
	)
	set( cxx
			-std=c++98
			-isystem external/boost_1_55_0/
			-isystem external/include/
			-isystem external/dbio/
	)
	set( compile
			-pipe
			-ffor-scope
	)
	set( warn
			-Wall
			-Wextra
			-pedantic
			-Werror
			-Wno-long-long
			-Wno-strict-aliasing
	)
endif()

# "gcc, 4.4"
if( ${COMPILER} STREQUAL "gcc" AND ${CMAKE_CXX_COMPILER_VERSION} MATCHES ".*4.4(.[0-9])*" )
	list( APPEND warn
			-Wno-uninitialized
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

# "gcc, 4.4, release"
if( ${COMPILER} STREQUAL "gcc" AND ${MODE} STREQUAL "release" AND ${CMAKE_CXX_COMPILER_VERSION} MATCHES ".*4.4(.[0-9])*" )
	list( REMOVE_ITEM compile
			-finline-limit=20000
	)
	list( APPEND compile
			-finline-limit=487
	)
endif()


###########################################################################
# Clang ###################################################################
###########################################################################

set(WITH_LIBSTDC++ ON CACHE BOOL "Build libraries using libstdc++ when using clang compiler. ")

if( ${COMPILER} STREQUAL "clang" )
	set( cc
			#-std=c99
			-isystem external/boost_1_55_0/
			-isystem external/include/
			-isystem external/dbio/
	)
	set( cxx
			-std=c++98
			-isystem external/boost_1_55_0/
			-isystem external/include/
			-isystem external/dbio/
	)
	set( compile
			-pipe
			-Qunused-arguments
			-DUNUSUAL_ALLOCATOR_DECLARATION
			-ftemplate-depth-256
	)
	set( warn
			-W
			-Wall
			-Wextra
			-pedantic
			#-Weverything
			-Werror
			-Wno-long-long
			-Wno-strict-aliasing
			#-Wno-documentation
			#-Wno-padded
			#-Wno-weak-vtables
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

    REMOVE_DEFINITIONS(-DPTR_MODERN)
    REMOVE_DEFINITIONS(-DPTR_BOOST)

    ADD_DEFINITIONS(-DCXX11)
    ADD_DEFINITIONS(-DPTR_STD)

    set( cxx
          -std=c++11
          -isystem external/boost_1_55_0/
          -isystem external/include/
          -isystem external/dbio/
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

# Pyrosetta build options
set(WITH_PYROSETTA OFF CACHE BOOL "Build rosetta libraries with pyrosetta-specific build configuration.")
if( WITH_PYROSETTA )
  list( APPEND defines -DPYROSETTA )
endif()

###########################################################################
# Extras ###################################################################
###########################################################################

if( EXTRAS )

	# cxx11"
	if( ${EXTRAS} STREQUAL "cxx11" )
	    set( cc
			    #-std=c99
			    -isystem external/boost_1_55_0/
			    -isystem external/include/
			    -isystem external/dbio/
	    )
	    set( cxx
			    -std=c++11
			    -isystem external/boost_1_55_0/
			    -isystem external/include/
			    -isystem external/dbio/

	    )
	    set( compile
			    -pipe
			    -ftemplate-depth-256
    			    -fPIC
			    -DBOOST_ERROR_CODE_HEADER_ONLY
			    -DBOOST_SYSTEM_NO_DEPRECATED
			    -I /usr/include
			    -I /usr/local/include
			    -I src
			    -I external/include
    	    )
	    if( ${COMPILER} STREQUAL "gcc" )
		list( APPEND compile
		    -I src/platform/linux
		)
	    endif()


	    if( ${COMPILER} STREQUAL "clang" )
		list( APPEND compile
		    -DUNUSUAL_ALLOCATOR_DECLARATION
		    -stdlib=libstdc++
		    -Qunused-arguments
		    -I src/platform/macos/64/clang/6.1
		    -I src/platform/macos/64/clang
		    -I src/platform/macos/64
		)

		set( shlink
		    -stdlib=libstdc++
		)

		set( link
		    -stdlib=libstdc++
		)

	    endif()

	endif()


	# "macos, graphics"
	if( APPLE AND ${EXTRAS} STREQUAL "graphics" )
		find_package( GLUT REQUIRED )
		find_package( OpenGL REQUIRED )
		include_directories( /usr/X11R6/include )
		link_directories( /usr/X11R6/lib )
		set( cxx
			    -std=c++98
			    -isystem external/boost_1_55_0/
			    -isystem external/include/
			    -isystem external/dbio/

	        )
		REMOVE_DEFINITIONS(-DCXX11)
		REMOVE_DEFINITIONS(-DPTR_STD)
		list( APPEND defines
				-DGL_GRAPHICS
				-DMAC
				-DPTR_MODERN
				-DPTR_BOOST
		)
		set( link
				-stdlib=libstdc++
				-framework GLUT
				-framework OpenGL
				-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib
		)
		set( shlink
				-stdlib=libstdc++
				-framework GLUT
				-framework OpenGL
				-dylib_file /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib
		)
		set( compile
				-pipe
				-ffast-math
				-funroll-loops
				-ftemplate-depth=250
				-stdlib=libstdc++
				-mmacosx-version-min=10.6
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
