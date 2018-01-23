CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(COMPILER clang)
SET(MODE release_bluegene)
SET(LINK_RULES static)

INCLUDE(../build/libraries.cmake)
INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/static.build.cmake)
INCLUDE(../build/build.settings.cmake)
#INCLUDE(../build/build.post.cmake)
INCLUDE(../build/external.static.cmake)

FOREACH(LIBRARY ${LIBRARIES})
	INCLUDE(../build/${LIBRARY}.cmake)
	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
	set_target_properties(${LIBRARY} PROPERTIES LINKER_LANGUAGE CXX)
	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
ENDFOREACH( LIBRARY )

LINK_DIRECTORIES(
	/bgsys/drivers/ppcfloor/comm/gcc/include/
	/soft/libraries/alcf/current/gcc/ZLIB/lib
)

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DUSEMPI)
ADD_DEFINITIONS(-DMULTI_THREADED)
ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
SET(COMPILE_FLAGS "-O3")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 " ) 

#FOREACH(LIBRARY ${LIBRARIES})
#	INCLUDE(../build/${LIBRARY}.cmake)
#	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
#	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
#		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
#	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
#ENDFOREACH( LIBRARY )
##### Libraries

### Static Libraries (including rosetta)
### These MUST be full paths, otherwise binaries will not recompile when these change
### THIS ORDER IS VERY IMPORTANT -> NEEDS TO BE REVERSE ORDER OF DEPENDENCY GRAPH
SET(LINK_STATIC_LIBS 
    ${CMAKE_CURRENT_SOURCE_DIR}/libdevel.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.8.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.7.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_e.6.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_d.6.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_c.6.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.6.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.6.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_h.5.a   
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_g.5.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_f.5.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_e.5.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_d.5.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_c.5.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.5.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.5.a           
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.4.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.3.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.2.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.2.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.1.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.5.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.4.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.3.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.2.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libcore.1.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libbasic.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libnumeric.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libutility.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libObjexxFCL.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libcppdb-static.a  
    ${CMAKE_CURRENT_SOURCE_DIR}/libcifparse.a  
    ${CMAKE_CURRENT_SOURCE_DIR}/liblibxml2.a  
    -lz) 
### libc must be linked dynamically if using static ompi libs
#SET(LINK_DYNAMIC_LIBS "-Wl,-Bdynamic -l c")

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic")

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic -lgcc_eh")
SET(TRAILING_STATIC "-Bstatic -L../../external/lib/")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ${TRAILING_STATIC})

### Change dynamic linker from default /lib/ld64.so.1 to /lib64/ld-2.5.so,
### since that file doesn't exist and path is hardcoded
SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--relax")



##### Top level binaries
# just loophash stuff for now
### loophash_createfiltereddb 
#ADD_EXECUTABLE(loophash_createfiltereddb.cxx11threadmpi.linuxbgclangrelease ../../src/apps/pilot/kenjung/loophash_createfiltereddb.cc)
#TARGET_LINK_LIBRARIES(loophash_createfiltereddb.cxx11threadmpi.linuxbgclangrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(loophash_createfiltereddb.cxx11threadmpi.linuxbgclangrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(loophash_createfiltereddb.cxx11threadmpi.linuxbgclangrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
#ADD_DEPENDENCIES(loophash_createfiltereddb.cxx11threadmpi.linuxbgclangrelease ObjexxFCL utility numeric basic 
#	core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols.6 protocols.7 devel)

### design_cycpeptide_MPI
#ADD_EXECUTABLE(design_cycpeptide_MPI.cxx11threadmpi.linuxbgclangrelease ../../src/apps/pilot/vmullig/design_cycpeptide_MPI.cc)
#TARGET_LINK_LIBRARIES(design_cycpeptide_MPI.cxx11threadmpi.linuxbgclangrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(design_cycpeptide_MPI.cxx11threadmpi.linuxbgclangrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(design_cycpeptide_MPI.cxx11threadmpi.linuxbgclangrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
#ADD_DEPENDENCIES(design_cycpeptide_MPI.cxx11threadmpi.linuxbgclangrelease cifparse libxml2 ObjexxFCL utility numeric
#basic core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2
#protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4
#protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4
#protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5
#protocols.6 protocols.7 devel)

### simple_cycpep_predict
ADD_EXECUTABLE(simple_cycpep_predict.cxx11threadmpi.linuxbgclangrelease
../../src/apps/public/cyclic_peptide/simple_cycpep_predict.cc)
TARGET_LINK_LIBRARIES(simple_cycpep_predict.cxx11threadmpi.linuxbgclangrelease "${LINK_ALL_LIBS}" )
SET_TARGET_PROPERTIES(simple_cycpep_predict.cxx11threadmpi.linuxbgclangrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(simple_cycpep_predict.cxx11threadmpi.linuxbgclangrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(simple_cycpep_predict.cxx11threadmpi.linuxbgclangrelease ObjexxFCL cifparse libxml2 utility numeric
basic core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2
protocols_b.2 protocols.3 protocols.4 protocols_a.5 protocols_b.5 protocols_c.5
protocols_d.5 protocols_e.5 protocols_f.5 protocols_g.5 protocols_h.5
protocols_a.6 protocols_b.6 protocols_c.6 protocols_d.6 protocols_e.6
protocols.7 protocols.8 devel)

### rosetta_scripts
ADD_EXECUTABLE(rosetta_scripts.cxx11threadmpi.linuxbgclangrelease ../../src/apps/public/rosetta_scripts/rosetta_scripts.cc)
TARGET_LINK_LIBRARIES(rosetta_scripts.cxx11threadmpi.linuxbgclangrelease "${LINK_ALL_LIBS}" )
SET_TARGET_PROPERTIES(rosetta_scripts.cxx11threadmpi.linuxbgclangrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(rosetta_scripts.cxx11threadmpi.linuxbgclangrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(rosetta_scripts.cxx11threadmpi.linuxbgclangrelease ObjexxFCL cifparse libxml2 utility numeric
basic core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2
protocols_b.2 protocols.3 protocols.4 protocols_a.5 protocols_b.5 protocols_c.5
protocols_d.5 protocols_e.5 protocols_f.5 protocols_g.5 protocols_h.5
protocols_a.6 protocols_b.6 protocols_c.6 protocols_d.6 protocols_e.6
protocols.7 protocols.8 devel)

### loophash_mpi 
#ADD_EXECUTABLE(loophash_mpi.cxx11threadmpi.linuxbgclangrelease ../../src/apps/pilot/mike/loophash_mpi.cc)
#TARGET_LINK_LIBRARIES(loophash_mpi.cxx11threadmpi.linuxbgclangrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(loophash_mpi.cxx11threadmpi.linuxbgclangrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(loophash_mpi.cxx11threadmpi.linuxbgclangrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
#ADD_DEPENDENCIES(loophash_mpi.cxx11threadmpi.linuxbgclangrelease ObjexxFCL utility numeric basic 
#	core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols.6 protocols.7 devel)

#### minirosetta 
#ADD_EXECUTABLE(minirosetta.cxx11threadmpi.linuxbgclangrelease ../../src/apps/public/boinc/minirosetta.cc)
#TARGET_LINK_LIBRARIES(minirosetta.cxx11threadmpi.linuxbgclangrelease ${LINK_ALL_LIBS})
#SET_TARGET_PROPERTIES(minirosetta.cxx11threadmpi.linuxbgclangrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(minirosetta.cxx11threadmpi.linuxbgclangrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
## Dependencies are set so that binaries are created after rosetta libraries
## are built
#ADD_DEPENDENCIES(minirosetta.cxx11threadmpi.linuxbgclangrelease ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols devel)
