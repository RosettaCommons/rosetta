OPTION(DEBUG "Build debugging symbols" OFF)
IF (DEBUG)
    SET(MODE debug)
ELSE ()
    SET(MODE release)
    ADD_DEFINITIONS(-DNDEBUG)
ENDIF ()
MESSAGE (STATUS "Build Type: ${MODE}")

