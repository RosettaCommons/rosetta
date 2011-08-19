# This mostly follows 
# URL: http://cppcms.svn.sourceforge.net/svnroot/cppcms/cppdb/trunk/CMakeLists.txt
# at revision -r1700

set(CPPDB_VERSION 0.1.0)
set(CPPDB_SOVERSION 0)


# General settings

include_directories(../../external/dbio)
include_directories(../../external)

option(DISABLE_SQLITE	"Link sqlite3 backend into the libcppdb" OFF)

add_definitions(-DCPPDB_EXPORTS)
add_definitions(-DCPPDB_LIBRARY_PREFIX="${CMAKE_SHARED_LIBRARY_PREFIX}")
add_definitions(-DCPPDB_LIBRARY_SUFFIX="${CMAKE_SHARED_LIBRARY_SUFFIX}")
add_definitions(-DCPPDB_SOVERSION="${CPPDB_SOVERSION}")
add_definitions(-DCPPDB_DISABLE_THREAD_SAFETY)
add_definitions(-DCPPDB_DISABLE_SHARED_OBJECT_LOADING)

# Backend configuration

set(INTERNAL_SOURCES)
set(INTERNAL_LIBRARIES)

if(NOT DISABLE_SQLITE)

	set(SQLITE3_SRC ../../external/dbio/sqlite3/sqlite3.c)
	add_definitions(-DSQLITE_OMIT_LOAD_EXTENSION)
	add_definitions(-DSQLITE_OMIT_DISABLE_LFS)
	add_definitions(-DSQLITE_THREADSAFE=0)

	add_library(sqlite3 SHARED ${SQLITE3_SRC})
        
        set(INTERNAL_SOURCES ${INTERNAL_SOURCES} ../../external/dbio/cppdb/sqlite3_backend.cpp)
        set(INTERNAL_LIBRARIES ${INTERNAL_LIBRARIES} sqlite3)
        add_definitions(-DCPPDB_WITH_SQLITE3)

endif()

# cppdb library configuration

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
	${INTERNAL_SOURCES}
	)

add_library(cppdb SHARED ${CPPDB_SRC})
set_target_properties(cppdb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}" )
set_target_properties(cppdb PROPERTIES LINK_FLAGS "${LINK_FLAGS}" )

foreach(LIB ${INTERNAL_LIBRARIES})
	target_link_libraries(cppdb ${LIB})
endforeach()

SET(LINK_EXTERNAL_LIBS ${LINK_EXTERNAL_LIBS} cppdb)