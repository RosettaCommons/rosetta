CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET( NACL_SDK_ROOT ../../external/native_client_sdk_0_2_803_0 )
SET( NACL_TOOLCHAIN_DIR toolchain/mac_x86  )
#SET( CMAKE_CXX_COMPILER  ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-g++ )
SET( CMAKE_LINKER        ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-g++ )
SET( CMAKE_AR            ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-ar ) 
SET( CMAKE_RANLIB        ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-ranlib )
SET( CMAKE_STRIP         ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-strip )

INCLUDE(../build/build.pre.cmake)
INCLUDE(../build/external.static.cmake)
INCLUDE( ${NACL_SDK_ROOT} )

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DNATCL)
SET(COMPILE_FLAGS " -Wall -Wno-long-long -pthread -Werror -O2")

SET( NACL_SEL_LDR32  ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl-sel_ldr )
SET( NACL_SEL_LDR64  ${NACL_SDK_ROOT}/${NACL_TOOLCHAIN_DIR}/bin/nacl64-sel_ldr )


SET( LIBRARIES zlib ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols )
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

#SET(LINK_STATIC_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.5.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.4.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.3.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.2.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.1.a ${CMAKE_CURRENT_SOURCE_DIR}/libbasic.a ${CMAKE_CURRENT_SOURCE_DIR}/libnumeric.a ${CMAKE_CURRENT_SOURCE_DIR}/libutility.a ${CMAKE_CURRENT_SOURCE_DIR}/libObjexxFCL.a ${CMAKE_CURRENT_SOURCE_DIR}/libcppdb-static.a   ${CMAKE_CURRENT_SOURCE_DIR}/libzlib.a -Xlinker --start-group -lppruntime -lppapi_cpp -lgoogle_nacl_platform -lgio -lpthread -lsrpc -lc -lg -liberty -lnosys  -Xlinker --end-group   -m32   ) 

SET(LINK_STATIC_LIBS   ${CMAKE_CURRENT_SOURCE_DIR}/libprotocols.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.5.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.4.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.3.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.2.a ${CMAKE_CURRENT_SOURCE_DIR}/libcore.1.a ${CMAKE_CURRENT_SOURCE_DIR}/libbasic.a ${CMAKE_CURRENT_SOURCE_DIR}/libnumeric.a ${CMAKE_CURRENT_SOURCE_DIR}/libutility.a ${CMAKE_CURRENT_SOURCE_DIR}/libObjexxFCL.a ${CMAKE_CURRENT_SOURCE_DIR}/libcppdb-static.a  ${CMAKE_CURRENT_SOURCE_DIR}/libzlib.a -Xlinker --start-group -lppruntime -lppapi_cpp -lgoogle_nacl_platform -lgio -lpthread -lsrpc -lc -lg -liberty -lnosys  -Xlinker --end-group   -m32   ) 


# trailing static needed so mpicc actually uses static ompi libs
#SET(TRAILING_STATIC "-Wl,-Bstatic")

### Put all the libraries together
### ORDER IS IMPORTANT, do not change unless you understand -Wl,-Bstatic
SET(LINK_ALL_LIBS ${LINK_STATIC_LIBS} ) # ${LINK_DYNAMIC_LIBS} ${TRAILING_STATIC})

### Change dynamic linker from default /lib/ld64.so.1 to /lib64/ld-2.5.so,
### since that file doesn't exist and path is hardcoded
SET(CMAKE_EXE_LINKER_FLAGS "  ")



##### Top level binaries
# just loophash stuff for now
### loophash_createfiltereddb 
ADD_EXECUTABLE(mini_rosetta_native_client.nexe ../../src/apps/pilot/mike/mini_rosetta_native_client.cc)
TARGET_LINK_LIBRARIES(mini_rosetta_native_client.nexe ${LINK_ALL_LIBS})
SET_TARGET_PROPERTIES(mini_rosetta_native_client.nexe PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(mini_rosetta_native_client.nexe PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
# Dependencies are set so that binaries are created after rosetta libraries
# are built
ADD_DEPENDENCIES(mini_rosetta_native_client.nexe zlib ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4 core.5 protocols )
#ADD_CUSTOM_COMMAND(
#    OUTPUT CMakeFiles/EXE.dir/link.txt
#		    COMMAND ${CMAKE_COMMAND} -E touch CMakeFiles/EXE.dir/link.txt
#				    DEPENDS linkerscript
#						)
#
