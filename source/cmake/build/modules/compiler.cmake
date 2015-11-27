IF (${CMAKE_CXX_COMPILER} MATCHES "clang")
    SET(COMPILER "clang")

    # Quiet warnings about non-clang compatable arguments (either gcc-specific ones, or newer-version ones.)
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} -Qunused-arguments")
ELSE ()
    SET(COMPILER "g++")
    
ENDIF ()

MESSAGE(STATUS "Compiler: ${COMPILER}")
