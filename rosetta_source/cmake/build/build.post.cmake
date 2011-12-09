#SET( LIBRARIES ObjexxFCL utility numeric basic core protocols devel )
SET( LIBRARIES
	ObjexxFCL utility numeric basic
	core.1 core.2 core.3 core.4 core.5 protocols.1 protocols_a.2 protocols_b.2 protocols.3 protocols_a.4 protocols_b.4 protocols_c.4 protocols_d.4 protocols_e.4 protocols_f.4 protocols_g.4 protocols_h.4 protocols_a.5 protocols_b.5 protocols_c.5 protocols_d.5 protocols_e.5 protocols_f.5 protocols.6 protocols.7 devel
)

## External built libraries
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
#FOREACH(LIBRARY ${LIBRARIES})
FOREACH(LIBRARY ${LIBRARIES})
	TARGET_LINK_LIBRARIES( ${LIBRARY} ${LINK_ROSETTA_LIBS} )
	SET( LINK_ROSETTA_LIBS "${LIBRARY};${LINK_ROSETTA_LIBS}" )
	ADD_DEPENDENCIES( BUILD_ROSETTA_LIBS ${LIBRARY} )
ENDFOREACH(LIBRARY)

SET(LINK_ALL_LIBS ${LINK_ROSETTA_LIBS} ${LINK_EXTERNAL_LIBS})
