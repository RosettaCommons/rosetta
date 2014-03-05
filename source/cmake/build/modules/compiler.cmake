IF (${CMAKE_CXX_COMPILER} MATCHES "clang")
    SET(COMPILER "clang")

    # Quiet warnings about non-clang compatable arguments (either gcc-specific ones, or newer-version ones.)
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Qunused-arguments")

    # Recent Macs can have issues with the default standard library - explicitly specify the GNU one.
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} -stdlib=libstdc++") 

ELSE ()
    SET(COMPILER "g++")
    
ENDIF ()

MESSAGE(STATUS "Compiler: ${COMPILER}")
