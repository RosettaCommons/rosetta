CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/external.static.cmake)

# ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-g)
#ADD_DEFINITIONS(-DUSEBOOSTSERIALIZE) # This is no longer used in Rosetta code
ADD_DEFINITIONS(-DUSEMPI) # only reason i enable this is for mpi rank on the tracer :)
ADD_DEFINITIONS(-DUSEBOOSTMPI)
ADD_DEFINITIONS(-DUSELUA)
ADD_DEFINITIONS("-std=c++0x")

# because we're using more libs
LINK_DIRECTORIES(../../external/lib/linux/amd64)

# and more includes
INCLUDE_DIRECTORIES(../../external/include/lua)
INCLUDE_DIRECTORIES(../../external/include/luabind)

SET(COMPILE_FLAGS "-pipe -w -O0")

SET( LIBRARIES
	ObjexxFCL utility numeric basic
	core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols.6 protocols.7 devel
    )
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
    -lz
    -lboost_serialization
    -lboost_mpi
    -lluabind
    -llua
    -ldl #libdl needed by liblua
    ) 
### libc must be linked dynamically if using static ompi libs
SET(LINK_DYNAMIC_LIBS "-Wl,-Bdynamic -lc -lpthread")

# trailing static needed so mpicc actually uses static ompi libs
SET(TRAILING_STATIC "-Wl,-Bstatic")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ${LINK_DYNAMIC_LIBS} ${TRAILING_STATIC})

### Change dynamic linker from default /lib/ld64.so.1 to /lib64/ld-2.5.so,
### since that file doesn't exist and path is hardcoded
SET(CMAKE_EXE_LINKER_FLAGS "-static -Wl,--dynamic-linker /lib64/ld-2.5.so")



##### Top level binaries
SET(BINARIES
    ../../src/apps/pilot/kenjung/pose_serialize_test.cc
    ../../src/apps/pilot/kenjung/lua_test.cc
    ../../src/apps/pilot/kenjung/elscripts.cc
)

FOREACH(BINARY ${BINARIES})
    GET_FILENAME_COMPONENT(BIN_NAME ${BINARY} NAME_WE)
    #SET(BIN_NAME ${BIN_NAME}.mpi)
    ADD_EXECUTABLE(${BIN_NAME} ${BINARY}) 
    TARGET_LINK_LIBRARIES(${BIN_NAME} ${LINK_ALL_LIBS})
    SET_TARGET_PROPERTIES(${BIN_NAME} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${BIN_NAME} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
    ADD_DEPENDENCIES(${BIN_NAME} ${LIBRARIES})
    #ADD_CUSTOM_COMMAND(TARGET ${BIN_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${BIN_NAME} /work/kenjung/minibin/)
ENDFOREACH( BINARY )
