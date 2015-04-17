INCLUDE(../build/libraries.cmake)
INCLUDE(../build/external.cmake)

## Library definitions
#SET( LAST_LIBRARY "empty" )
FOREACH(LIBRARY ${LIBRARIES})
	INCLUDE(../build/${LIBRARY}.cmake)
	ADD_LIBRARY(${LIBRARY} ${LINK_TYPE} ${${LIBRARY}_files})
	IF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
		ADD_DEPENDENCIES( ${project} ${LAST_LIBRARY} )
	ENDIF( ${LAST_LIBRARY} NOT STREQUAL "empty" )
ENDFOREACH( LIBRARY )

## Library link order
ADD_CUSTOM_TARGET( BUILD_ROSETTA_LIBS )
SET( LINK_ROSETTA_LIBS ${LINK_EXTERNAL_LIBS} )
SET( LINK_LAST_LIB ${LINK_EXTERNAL_LIBS} )
FOREACH(LIBRARY ${LIBRARIES})
	TARGET_LINK_LIBRARIES( ${LIBRARY} ${LINK_LAST_LIB} )
	SET( LINK_ROSETTA_LIBS "${LIBRARY};${LINK_ROSETTA_LIBS}" )
	ADD_DEPENDENCIES( BUILD_ROSETTA_LIBS ${LIBRARY} )
	SET( LINK_LAST_LIB ${LIBRARY} )
ENDFOREACH(LIBRARY)

# The following links the external libs twice, which is weird even if it doesn't cause
# problems.
#SET(LINK_ALL_LIBS "${LINK_ROSETTA_LIBS};${LINK_EXTERNAL_LIBS}")
SET(LINK_ALL_LIBS ${LINK_ROSETTA_LIBS})

# Add custom targets for just the applications (Dependancies added elsewhere)
ADD_CUSTOM_TARGET(apps DEPENDS BUILD_ROSETTA_LIBS) # Just released applications
ADD_CUSTOM_TARGET(pilot_apps DEPENDS BUILD_ROSETTA_LIBS) # Just pilot applications
ADD_CUSTOM_TARGET(bin DEPENDS apps pilot_apps) #All executables
