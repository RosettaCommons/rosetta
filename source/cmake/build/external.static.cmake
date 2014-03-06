# This mostly follows 
# URL: http://cppcms.svn.sourceforge.net/svnroot/cppcms/cppdb/trunk/CMakeLists.txt
# at revision -r1700

set(CPPDB_VERSION 0.1.0)
set(CPPDB_SOVERSION 0)


# General settings

include_directories(SYSTEM ../../external/dbio/sqlite3)
include_directories(SYSTEM ../../external/dbio)
include_directories(SYSTEM ../../external)

option(DISABLE_SQLITE	"Link sqlite3 backend into the libcppdb" OFF)

add_definitions(-DCPPDB_EXPORTS)
add_definitions(-DCPPDB_LIBRARY_PREFIX="${CMAKE_SHARED_LIBRARY_PREFIX}")
add_definitions(-DCPPDB_LIBRARY_SUFFIX="${CMAKE_SHARED_LIBRARY_SUFFIX}")
add_definitions(-DCPPDB_SOVERSION="${CPPDB_SOVERSION}")
add_definitions(-DCPPDB_DISABLE_THREAD_SAFETY)
add_definitions(-DCPPDB_DISABLE_SHARED_OBJECT_LOADING)
add_definitions(-DSQLITE_OMIT_LOAD_EXTENSION)
add_definitions(-DSQLITE_OMIT_DISABLE_LFS)
add_definitions(-DSQLITE_THREADSAFE=0)
add_definitions(-DCPPDB_MAJOR=0)
add_definitions(-DCPPDB_MINOR=3)
add_definitions(-DCPPDB_PATCH=0)
add_definitions(-DCPPDB_VERSION="0.3.0")

# Backend configuration

set(INTERNAL_SOURCES)
set(INTERNAL_LIBRARIES)

if(NOT DISABLE_SQLITE)

	# This is currenlty using the system version of sqlite3.
        # Rosetta distributes sqlite3 source files in
        # external/dbio/sqlite3, so the scons build system can build
        # it from source for consistency.  Most modern systems,
        # however, have sqlite3 already, so not having to build it
        # will be slightly faster.  If you are on a system without
        # sqlite3, you may want to build the one from source.  See
        # external/SConscript.external to see how it is built via the
        # scons build system.

	find_library(SQLITE3_LIB sqlite3)
	find_path(SQLITE3_PATH sqlite3.h)
	if(SQLITE3_LIB AND SQLITE3_PATH)
		include_directories(SYSTEM ${SQLITE3_PATH})
		add_definitions(-DCPPDB_WITH_SQLITE3)
		set(INTERNAL_SOURCES ${INTERNAL_SOURCES} ../../external/dbio/cppdb/sqlite3_backend.cpp)
		set(INTERNAL_LIBRARIES ${INTERNAL_LIBRARIES} ${SQLITE3_LIB})
	else()
		message("-- sqlite3 library was not found, disabling sqlite3 backend")
	endif()
endif()

# Main library configuration

set(CPPDB_SRC
	../../external/dbio/cppdb/utils.cpp
	../../external/dbio/cppdb/mutex.cpp
	../../external/dbio/cppdb/driver_manager.cpp
	../../external/dbio/cppdb/conn_manager.cpp
	../../external/dbio/cppdb/shared_object.cpp
	../../external/dbio/cppdb/pool.cpp
	../../external/dbio/cppdb/backend.cpp
	../../external/dbio/cppdb/frontend.cpp
	../../external/dbio/cppdb/atomic_counter.cpp
	../../external/dbio/sqlite3/sqlite3.c 
	${INTERNAL_SOURCES}
	)

#add_library(cppdb SHARED ${CPPDB_SRC})
add_library(cppdb-static STATIC ${CPPDB_SRC})
set_target_properties(cppdb-static PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}" )
set_target_properties(cppdb-static PROPERTIES LINK_FLAGS "${LINK_FLAGS}" )

foreach(LIB ${INTERNAL_LIBRARIES})
	target_link_libraries(cppdb-static ${LIB})
endforeach()

SET(LINK_EXTERNAL_LIBS ${LINK_EXTERNAL_LIBS} cppdb-static)
