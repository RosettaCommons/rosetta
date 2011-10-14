CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/external.static.cmake)

LINK_DIRECTORIES(/bgsys/drivers/ppcfloor/arch/include /soft/apps/zlib-1.2.3/lib)

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DUSEMPI)
ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
SET(COMPILE_FLAGS "-O3")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 " ) 

SET( LIBRARIES ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols devel )
FOREACH(LIBRARY ${LIBRARIES})
	INCLUDE(../build/${LIBRARY}.cmake)
	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
ENDFOREACH( LIBRARY )
##### Libraries

### Static Libraries (including rosetta)
### These MUST be full paths, otherwise binaries will not recompile when these change
### THIS ORDER IS VERY IMPORTANT -> NEEDS TO BE REVERSE ORDER OF DEPENDENCY GRAPH
SET(LINK_STATIC_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/libdevel.a ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.5.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.4.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.3.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.2.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.1.a ${CMAKE_CURRENT_SOURCE_DIR}/libbasic.a ${CMAKE_CURRENT_SOURCE_DIR}/libnumeric.a ${CMAKE_CURRENT_SOURCE_DIR}/libutility.a ${CMAKE_CURRENT_SOURCE_DIR}/libObjexxFCL.a ${CMAKE_CURRENT_SOURCE_DIR}/libcppdb-static.a  -lz) 
### libc must be linked dynamically if using static ompi libs
#SET(LINK_DYNAMIC_LIBS "-Wl,-Bdynamic -l c")

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
#SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ${LINK_DYNAMIC_LIBS} ${TRAILING_STATIC})

SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ) 

### Change dynamic linker from default /lib/ld64.so.1 to /lib64/ld-2.5.so,
### since that file doesn't exist and path is hardcoded
SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl")



##### Top level binaries
# just loophash stuff for now
### loophash_createfiltereddb 
ADD_EXECUTABLE(loophash_createfiltereddb.mpi.linuxgccrelease ../../src/apps/pilot/kenjung/loophash_createfiltereddb.cc)
TARGET_LINK_LIBRARIES(loophash_createfiltereddb.mpi.linuxgccrelease ${LINK_ALL_LIBS})
SET_TARGET_PROPERTIES(loophash_createfiltereddb.mpi.linuxgccrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(loophash_createfiltereddb.mpi.linuxgccrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(loophash_createfiltereddb.mpi.linuxgccrelease ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols devel)

### loophash_mpi 
ADD_EXECUTABLE(loophash_mpi.mpi.linuxgccrelease ../../src/apps/pilot/mike/loophash_mpi.cc)
TARGET_LINK_LIBRARIES(loophash_mpi.mpi.linuxgccrelease ${LINK_ALL_LIBS})
SET_TARGET_PROPERTIES(loophash_mpi.mpi.linuxgccrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(loophash_mpi.mpi.linuxgccrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(loophash_mpi.mpi.linuxgccrelease ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols devel)

### minirosetta 
ADD_EXECUTABLE(minirosetta.mpi.linuxgccrelease ../../src/apps/public/boinc/minirosetta.cc)
TARGET_LINK_LIBRARIES(minirosetta.mpi.linuxgccrelease ${LINK_ALL_LIBS})
SET_TARGET_PROPERTIES(minirosetta.mpi.linuxgccrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(minirosetta.mpi.linuxgccrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(minirosetta.mpi.linuxgccrelease ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols devel)
