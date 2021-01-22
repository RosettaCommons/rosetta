# The settings below are combined together to create the final set of settings
# for CMake + Ninja builds. These settings attempt to match those found in
# tools/build/basic.settings, which are used for building with scons.
#
# Optional extra dependencies, which may be activated via cmake options, are
# included in this file. This may correspond to extras defined in
# basic.settings.
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

	if( ${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "8.1" )
		list( APPEND compile
			-Wl,--no-as-needed
		)
		list( APPEND link
			-Wl,--no-as-needed
		)
		list( REMOVE_ITEM compile
			-ffor-scope
		)
	endif()

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
			-D_GLIBCXX_DEBUG
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

# "gcc, release_static"
if( ${COMPILER} STREQUAL "gcc" AND ${MODE} STREQUAL "release_static" )
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
	list( REMOVE_ITEM warn
			-Werror
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
			-I ../../external/boost_submod
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
	#list( APPEND mode
	#		-g
	#		-ggdb
	#)
endif()

###########################################################################
# Windows Visual Studio - CL ##############################################
###########################################################################

if( ${COMPILER} STREQUAL "cl" )
	ADD_DEFINITIONS(-DPTR_STD)
	ADD_DEFINITIONS(-DCXX11)
	set( cc
			# MSVC doesn't support setting the C version to use.
	)
	set( cxx
			# MSVC uses C++14 as the default, there is no way to request C++11.
			#-std=c++11
	)
	set( compile
		-I../../src/platform/windows/msvc/

		# Use the updated MSVC frontend that has improved standards conformance.
		-permissive-
		# Enable >64k sections in obj files.
		-bigobj

		# Disable noisy warnings.
		# Signed/unsigned mismatch for comparison.
		-wd4018
		# Unreferenced local variable.
		-wd4101
		# "*/" found outside of comment.
		-wd4138
		# Conversion from floating point to integer may lose data.
		-wd4244
		# Ensure class has a dll-interface if exported.
		-wd4251
		# Conversion from size_t to smaller type may lose data.
		-wd4267
		# Non dll-exported class used as the base for dll-exported class.
		-wd4275
		# Conversion to smaller type during initialization/construction may lose data.
		-wd4305
		# Use of insecure CRT function.
		-wd4996
	)

	# Prevent the Windows headers from defining the "min" and "max" macros.
	ADD_DEFINITIONS(-DNOMINMAX)
	# Use the reduced version of the Windows headers (e.g., no Winsock 1, no OLE2, etc.).
	ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)

	# Force Boost.System to be header-only.
	ADD_DEFINITIONS(-DBOOST_ERROR_CODE_HEADER_ONLY)
	# Remove deprecated features in Boost.System.
	ADD_DEFINITIONS(-DBOOST_SYSTEM_NO_DEPRECATED)
	# Disable automatic linking for Boost.
	ADD_DEFINITIONS(-DBOOST_ALL_NO_LIB)

	# Required to build rdkit on MSVC (see source\external\rdkit\GraphMol\SLNParse\CMakeLists.txt)
	ADD_DEFINITIONS(-DYY_NO_UNISTD_H)

	# Only permit four link commands to be run at once.
	# Rosetta builds *massive* static libraries which are linked together to make executables,
	# these libraries are massive enough that they can cause MSVC's linker comsume all physical
	# memory if too instances of link.exe are running at once, so prefer build reliability over
	# speed by forcing all link commands into a ninja job pool with max parallelism=4.
	# Additionally, link.exe is multi-threaded so it is best not to create as many link processes
	# as there are cores on the machine.
	set_property(GLOBAL PROPERTY JOB_POOLS link_job_pool=4)
	set(CMAKE_JOB_POOL_LINK link_job_pool)
endif()


# modes ###################################################################

# "cl, debug"
if( ${COMPILER} STREQUAL "cl" AND ${MODE} STREQUAL "debug" )
	list( APPEND compile
		# Place functions in individual sections (allows functions to be dropped by opt:ref).
		-Gy
	)
	list( APPEND shlink
		# Enable dropping unused functions/data - this is disabled by default with /DEBUG.
		-opt:ref
	)
	list( APPEND exelink
		# Enable dropping unused functions/data - this is disabled by default with /DEBUG.
		-opt:ref
	)

	foreach(t EXE SHARED)
		# Disable incremental linking (the programs linked in Rosetta are too large to be incrementally linked).
		string(REPLACE "/INCREMENTAL" "/INCREMENTAL:NO" CMAKE_${t}_LINKER_FLAGS_DEBUG "${CMAKE_${t}_LINKER_FLAGS_DEBUG}")

		# Enable fastlink for debug symbols (.pdb files) - this causes the PDB files for libraries and executables
		# to point back to the PDB files for .objs or .libs instead of copying their contents: greatly reducing the
		# time it takes to do linking and the size of the generated PDB files.
		string(REPLACE "/debug" "/debug:fastlink" CMAKE_${t}_LINKER_FLAGS_DEBUG "${CMAKE_${t}_LINKER_FLAGS_DEBUG}")
	endforeach()
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
	list( APPEND defines

			# disbaling due to linking problem with Xcode 9.0 -D_LIBCPP_DEBUG
			-D_GLIBCXX_DEBUG
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
			-I ../../external/boost_submod
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
	#list( APPEND mode
	#		-g
	#		-ggdb
	#)
endif()

###########################################################################
# Optional build components
###########################################################################

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../build/finders)
INCLUDE(../build/modules/mysql.cmake)
INCLUDE(../build/modules/hdf5.cmake)

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


	# serialization build
	if( ${EXTRAS} MATCHES "serialization" )
	    ADD_DEFINITIONS(-DSERIALIZATION)

	    # list( REMOVE_ITEM warn
	    # 	 -Werror
    	    # 	 -Wunused-function
	    # )

	    list ( APPEND warn
	    	# -Wno-deprecated
	    	# -Wno-unused-variable
	    	# -Wno-unused-parameter
	    	# -Wno-type-limits
    		-Wno-unused-function
	    )


	endif()

	if( ${EXTRAS} MATCHES "zeromq" )
	    ADD_DEFINITIONS(-DSERIALIZATION -DZEROMQ)

	    list ( APPEND warn
    		-Wno-unused-function
	    )

	    list( APPEND compile
			-pthread
	    )

	    list( APPEND link
			-lpthread
	    )

	    list( APPEND shlink
			-lpthread
	    )

	endif()

	if( ${EXTRAS} MATCHES "cxx11thread" )
	    ADD_DEFINITIONS(-DMULTI_THREADED)

	    list( APPEND compile
			-pthread
	    )

	    list( APPEND link
			-lpthread
	    )

	    list( APPEND shlink
			-lpthread
	    )
	endif()


endif()

# Make sure that the submodules are up-to-date w/r/t the extras
if( WIN32 )
	# EXECUTE_PROCESS on Windows requires that the "child process" is actually a binary executable, hence we use cmd to then invoke our
	# script. Additionally, cmd can only handle back slashes...
	EXECUTE_PROCESS(COMMAND "cmd" "/c ..\\..\\update_submodules.sh ${EXTRAS}" RESULT_VARIABLE rv )
else()
	EXECUTE_PROCESS(COMMAND "../../update_submodules.sh" "${EXTRAS}" RESULT_VARIABLE rv)
endif()
if( NOT rv STREQUAL "0" )
	message( FATAL_ERROR "update_submodules.sh failed with: ${rv}" )
endif()


###########################################################################
###########################################################################
###########################################################################

add_definitions( ${defines} )

foreach( flag ${shlink} )
	set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${flag}" )
endforeach()

if( WIN32 )
	foreach( flag ${exelink} )
		set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${flag}" )
	endforeach()
endif()

foreach( flag ${cc} ${cxx} ${compile} ${warn} ${mode} )
	set( COMPILE_FLAGS "${COMPILE_FLAGS} ${flag}" )
endforeach()

foreach( flag ${COMPILE_FLAGS} )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" )
endforeach()
