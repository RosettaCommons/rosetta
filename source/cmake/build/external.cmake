INCLUDE(../build/external_libraries.cmake)

## Library definitions
FOREACH(LIBRARY ${EXTERNAL_LIBRARIES})
	INCLUDE(../build/external_${LIBRARY}.cmake)
	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
	SET_TARGET_PROPERTIES(${LIBRARY} PROPERTIES COMPILE_DEFINITIONS "${${LIBRARY}_defines}")
	# Because CMAKE can't actually use CMAKE-formated lists...
	STRING( REPLACE ";" " " linkflags_string "${${LIBRARY}_linkflags}" )
	SET_TARGET_PROPERTIES(${LIBRARY} PROPERTIES LINK_FLAGS "${linkflags_string}")
	# Because CMAKE can't actually use CMAKE-formated lists...
	STRING( REPLACE ";" " " compileflags_string "${${LIBRARY}_compileflags}" )
	SET_TARGET_PROPERTIES(${LIBRARY} PROPERTIES COMPILE_FLAGS "${compileflags_string} -Wno-error -w") #Also Turn off warnings for external libraries
	SET(LINK_EXTERNAL_LIBS ${LINK_EXTERNAL_LIBS} ${LIBRARY})
ENDFOREACH( LIBRARY )

# This mostly follows
# URL: http://cppcms.svn.sourceforge.net/svnroot/cppcms/cppdb/trunk/CMakeLists.txt
# at revision -r1700

set(CPPDB_SOVERSION 0)

# General settings
# We use the SYSTEM directive here to avoid printing warning about code in the external headers
# (Warnings in the .cc and .cpp files themselves will still be printed.)
include_directories(SYSTEM ../../external/dbio/sqlite3)
include_directories(SYSTEM ../../external/dbio)
include_directories(SYSTEM ../../external)
include_directories(SYSTEM ../../external/libxml2/include)

option(DISABLE_SQLITE	"Link sqlite3 backend into the libcppdb" OFF)

add_definitions(-DCPPDB_EXPORTS)
add_definitions(-DCPPDB_LIBRARY_PREFIX="${CMAKE_SHARED_LIBRARY_PREFIX}")
add_definitions(-DCPPDB_LIBRARY_SUFFIX="${CMAKE_SHARED_LIBRARY_SUFFIX}")
add_definitions(-DCPPDB_SOVERSION="${CPPDB_SOVERSION}")
add_definitions(-DCPPDB_MAJOR=0)
add_definitions(-DCPPDB_MINOR=3)
add_definitions(-DCPPDB_PATCH=0)
add_definitions(-DCPPDB_VERSION="0.3.0")

# I don't know why thread safety and shared object loading were disabled, but I
# had to reenable both options in order to use mysql from an MPI job.
# kalekundert@ucsf.edu 2/10/14
#add_definitions(-DCPPDB_DISABLE_THREAD_SAFETY)
#add_definitions(-DCPPDB_DISABLE_SHARED_OBJECT_LOADING)

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

# set(LIBXML2_SRC
# 	../../external/libxml2/buf.c
# 	../../external/libxml2/c14n.c
# 	../../external/libxml2/catalog.c
# 	../../external/libxml2/chvalid.c
# 	../../external/libxml2/debugXML.c
# 	../../external/libxml2/dict.c
# 	../../external/libxml2/DOCBparser.c
# 	../../external/libxml2/encoding.c
# 	../../external/libxml2/entities.c
# 	../../external/libxml2/error.c
# 	../../external/libxml2/globals.c
# 	../../external/libxml2/hash.c
# 	../../external/libxml2/HTMLparser.c
# 	../../external/libxml2/HTMLtree.c
# 	../../external/libxml2/legacy.c
# 	../../external/libxml2/list.c
# 	../../external/libxml2/nanoftp.c
# 	../../external/libxml2/nanohttp.c
# 	../../external/libxml2/parser.c
# 	../../external/libxml2/parserInternals.c
# 	../../external/libxml2/pattern.c
# 	../../external/libxml2/relaxng.c
# 	../../external/libxml2/SAX2.c
# 	../../external/libxml2/SAX.c
# 	../../external/libxml2/schematron.c
# 	../../external/libxml2/threads.c
# 	../../external/libxml2/tree.c
# 	../../external/libxml2/uri.c
# 	../../external/libxml2/valid.c
# 	../../external/libxml2/xinclude.c
# 	../../external/libxml2/xlink.c
# 	../../external/libxml2/xmlIO.c
# 	../../external/libxml2/xmlmemory.c
# 	../../external/libxml2/xmlmodule.c
# 	../../external/libxml2/xmlreader.c
# 	../../external/libxml2/xmlregexp.c
# 	../../external/libxml2/xmlsave.c
# 	../../external/libxml2/xmlschemas.c
# 	../../external/libxml2/xmlschemastypes.c
# 	../../external/libxml2/xmlstring.c
# 	../../external/libxml2/xmlunicode.c
# 	../../external/libxml2/xmlwriter.c
# 	../../external/libxml2/xpath.c
# 	../../external/libxml2/xpointer.c
# 	../../external/libxml2/xzlib.c
# 	)
# 
# 	add_library( libxml2 SHARED ${LIBXML2_SRC} )


# "gcc, 4.9"
if( ${COMPILER} STREQUAL "gcc" AND ${CMAKE_CXX_COMPILER_VERSION} MATCHES ".*4[.]9[.]([0-9])*" )
  # Fix cppdb linking error in recent versions of g++
  set_target_properties(cppdb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -fno-strict-aliasing -ldl" )
  set_target_properties(cppdb PROPERTIES LINK_FLAGS "${LINK_FLAGS} -Wl,--no-as-needed -ldl" )
# "gcc, 5.x, 6.x"
elseif( ${COMPILER} STREQUAL "gcc" AND ${CMAKE_CXX_COMPILER_VERSION} MATCHES ".*[56][.]([0-9])*[.]([0-9])*" )
  set_target_properties(cppdb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -fno-strict-aliasing -ldl" )
  set_target_properties(cppdb PROPERTIES LINK_FLAGS "${LINK_FLAGS} -Wl,--no-as-needed -ldl" )
else()
  set_target_properties(cppdb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}" )
  set_target_properties(cppdb PROPERTIES LINK_FLAGS "${LINK_FLAGS}" )
endif()

foreach(LIB ${INTERNAL_LIBRARIES})
	target_link_libraries(cppdb ${LIB})
	#target_link_libraries(libxml2 ${LIB})
endforeach()

SET(LINK_EXTERNAL_LIBS ${LINK_EXTERNAL_LIBS} cppdb)
#SET(LINK_EXTERNAL_LIBS ${LINK_EXTERNAL_LIBS} libxml2)


