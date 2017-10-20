CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(COMPILER gcc)
SET(MODE release_static)
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
	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
ENDFOREACH( LIBRARY )

ADD_DEFINITIONS(-DNDEBUG)

### Static Libraries (including rosetta)
### These MUST be full paths, otherwise binaries will not recompile when these change
### THIS ORDER IS VERY IMPORTANT -> NEEDS TO BE REVERSE ORDER OF DEPENDENCY GRAPH
SET(LINK_STATIC_LIBS 
    ${CMAKE_CURRENT_SOURCE_DIR}/libdevel.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.7.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.6.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_e.5.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_d.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_c.5.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.5.a
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.5.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_h.4.a   
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_g.4.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_f.4.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_e.4.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_d.4.a 
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_c.4.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_b.4.a    
    ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols_a.4.a           
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

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic -lgcc_eh")
SET(TRAILING_STATIC "-Bstatic -L../../external/lib/")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ${TRAILING_STATIC})

### since that file doesn't exist and path is hardcoded
SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--relax")


### PUT THE BINARIES TO BUILD HERE!
SET( BINARIES
	../../src/apps/public/rosetta_scripts/rosetta_scripts.cc
	../../src/apps/public/relax.cc
	../../src/apps/public/analysis/score_jd2.cc
	../../src/apps/public/extract_pdbs.cc
)

# Some intense makefile string parsing here... Hurt my soul to learn this.
FOREACH(BINARY ${BINARIES})
	string(REPLACE "/" ";" BINARY_LIST ${BINARY})
	list(GET BINARY_LIST -1 BASE_FN)
	string(REPLACE "." ";" BASE_FN_LIST ${BASE_FN})
	list(GET BASE_FN_LIST 0 FN_NO_EXT)
	SET(BINARY_NAME ${FN_NO_EXT}.static.linux${COMPILER}release)
	ADD_EXECUTABLE(${BINARY_NAME} ${BINARY})
	TARGET_LINK_LIBRARIES(${BINARY_NAME} "${LINK_ALL_LIBS}" )
	SET_TARGET_PROPERTIES(${BINARY_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
	SET_TARGET_PROPERTIES(${BINARY_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
	ADD_DEPENDENCIES(${BINARY_NAME} ObjexxFCL cifparse libxml2 utility numeric
		basic core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2
		protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4
		protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4
		protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5
		protocols.6 protocols.7 devel
	)
ENDFOREACH( BINARY )

# TEMPLATE FOR ^^
#### rosetta_scripts
# ADD_EXECUTABLE(rosetta_scripts.static.linuxgccrelease ../../src/apps/public/rosetta_scripts/rosetta_scripts.cc)
# TARGET_LINK_LIBRARIES(rosetta_scripts.static.linuxgccrelease "${LINK_ALL_LIBS}" )
# SET_TARGET_PROPERTIES(rosetta_scripts.static.linuxgccrelease PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
# SET_TARGET_PROPERTIES(rosetta_scripts.static.linuxgccrelease PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# # Dependencies are set so that binaries are created after rosetta libraries
# # are built
# ADD_DEPENDENCIES(rosetta_scripts.static.linuxgccrelease ObjexxFCL cifparse libxml2 utility numeric
# basic core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2
# protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4
# protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4
# protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5
# protocols.6 protocols.7 devel)
