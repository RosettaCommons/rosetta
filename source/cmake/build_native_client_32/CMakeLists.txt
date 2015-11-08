CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#/Users/mtyka/src/nacl_sdk/pepper_23/toolchain/mac_x86_newlib/bin/i686-nacl-g++
SET( NACL_SDK_ROOT /Users/mtyka/src/nacl_sdk/pepper_23/ ) #../../external/native_client_sdk_0_2_803_0 )
SET( NACL_TOOLCHAIN_DIR toolchain/mac_x86_newlib  )
#SET( CMAKE_CC_COMPILER   ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/i686-nacl-gcc )
#SET( CMAKE_CXX_COMPILER  ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/i686-nacl-g++ )
SET( CMAKE_LINKER        ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/i686-nacl-g++ )
SET( CMAKE_AR            ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/i686-nacl-ar ) 
SET( CMAKE_RANLIB        ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/i686-nacl-ranlib )
SET( CMAKE_STRIP         ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/i686-nacl-strip )

SET( DISABLE_SQLITE )

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/external.static.cmake)
INCLUDE_DIRECTORIES( ${NACL_SDK_ROOT}/include/ )

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-D__native_client__)
ADD_DEFINITIONS(-DUSE_FILE_PROVIDER)
ADD_DEFINITIONS(-DEXIT_THROWS_EXCEPTION)
ADD_DEFINITIONS(-O2)
ADD_DEFINITIONS(-s)
ADD_DEFINITIONS(-m32)
ADD_DEFINITIONS(-MMD)
ADD_DEFINITIONS(-pthread)
ADD_DEFINITIONS(-Wno-long-long)
ADD_DEFINITIONS(-Wall)
ADD_DEFINITIONS(-Wswitch-enum)
ADD_DEFINITIONS(-pedantic)



#SET( NACL_SEL_LDR32  ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-sel_ldr )
#SET( NACL_SEL_LDR64  ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl64-sel_ldr )


SET( LIBRARIES
  ObjexxFCL
  utility
  numeric
  basic
  core.1
  core.2
  core.3
  core.4
  core.5
  protocols.1
  protocols_a.2
  protocols_b.2
  protocols.3
  protocols_a.4
  protocols_b.4
  protocols_c.4
  protocols_d.4
  protocols_e.4
  protocols_f.4
  protocols_g.4
  protocols_h.4
  protocols_a.5
  protocols_b.5
  protocols_c.5
  protocols_d.5
  protocols_e.5
  protocols.6
  protocols.7
#  devel
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
#    ${CMAKE_CURRENT_SOURCE_DIR}/libdevel.a
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
    ${CMAKE_CURRENT_SOURCE_DIR}/libz-static.a  
    --start-group -lnosys -lppapi -lppapi_cpp -lplatform -lgio -lpthread -lc -lcrt_common   --end-group ) 

# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ) # ${LINK_DYNAMIC_LIBS} ${TRAILING_STATIC})

### Change dynamic linker from default /lib/ld64.so.1 to /lib64/ld-2.5.so,
### since that file doesn't exist and path is hardcoded

SET(CMAKE_C_LINK_FLAGS "-Wl,-as-needed  -m32 -pthread " )
SET(CMAKE_CXX_LINK_FLAGS "-Wl,-as-needed  -m32 -pthread " )


##### Top level binaries
# just loophash stuff for now
### loophash_createfiltereddb 
ADD_EXECUTABLE(rosetta_native_client.nexe ../../src/apps/pilot/mike/rosetta_native_client.cc)
TARGET_LINK_LIBRARIES(rosetta_native_client.nexe ${LINK_ALL_LIBS})
SET_TARGET_PROPERTIES(rosetta_native_client.nexe PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(rosetta_native_client.nexe PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(rosetta_native_client.nexe zlib ${LIBRARIES} ) 
#ADD_CUSTOM_COMMAND(
#    OUTPUT CMakeFiles/EXE.dir/link.txt
#		    COMMAND ${CMAKE_COMMAND} -E touch CMakeFiles/EXE.dir/link.txt
#				    DEPENDS linkerscript
#						)
#



