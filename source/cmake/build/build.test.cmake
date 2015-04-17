INCLUDE(../build/test_libraries.cmake)

INCLUDE_DIRECTORIES(../../test)

ADD_CUSTOM_TARGET(unit DEPENDS BUILD_ROSETTA_LIBS)

## Library definitions
FOREACH(LIBRARY ${TEST_LIBRARIES})
	INCLUDE(../build/test_${LIBRARY}.cmake)
	#Need directories present for cxxtestgen.py to work
	FOREACH(DIRECTORY ${${LIBRARY}_testdirectories} )
		FILE( MAKE_DIRECTORY ${DIRECTORY} )
	ENDFOREACH( DIRECTORY )

	#Copy over input files, but only if changed
	FOREACH(INPUTFILE ${${LIBRARY}_testinputs} )
		ADD_CUSTOM_COMMAND(OUTPUT ${INPUTFILE} COMMAND ${CMAKE_COMMAND} -E copy ${SRCDIR}${INPUTFILE} ${INPUTFILE} MAIN_DEPENDENCY ${SRCDIR}${INPUTFILE} )
	ENDFOREACH( INPUTFILE )
	ADD_CUSTOM_TARGET( ${LIBRARY}.testinputs ALL DEPENDS ${${LIBRARY}_testinputs} )

	#Generate sub-test files
	FOREACH(TESTFILE ${${LIBRARY}_testfiles})
		STRING(REGEX REPLACE "[.]cc" ".hh" SOURCE ${SRCDIR}${TESTFILE} )
		SET( TARGET ${TESTFILE} )
		GET_FILENAME_COMPONENT( DIRECTORY ${TARGET} PATH )
		ADD_CUSTOM_COMMAND(OUTPUT ${TARGET} COMMAND ../../external/cxxtest/cxxtestgen.py --have-std --part -o ${TARGET} ${SOURCE} MAIN_DEPENDENCY ${SOURCE} )
 	ENDFOREACH( TESTFILE )

	#Generate main wrapper application for library test
	SET( TARGET ${LIBRARY}.cxxtest.cpp )
	SET( SOURCE ${SRCDIR}${LIBRARY}/${LIBRARY}.cxxtest.hh )
	ADD_CUSTOM_COMMAND(OUTPUT ${TARGET} COMMAND ../../external/cxxtest/cxxtestgen.py --error-printer --root -o ${TARGET} ${SOURCE} MAIN_DEPENDENCY ${SOURCE} )
	ADD_EXECUTABLE( ${LIBRARY}.test ${LIBRARY}.cxxtest.cpp ${${LIBRARY}_testfiles} )
	# There's probably a smarter way of just linking in those libraries which are relevant,
	# but with splitting out core & protocols into many sub-libraries, I'm not sure what it is.
	TARGET_LINK_LIBRARIES( ${LIBRARY}.test ${LINK_ALL_LIBS} ) #  ${${LIBRARY}_testfiles}

	ADD_DEPENDENCIES( unit ${LIBRARY}.test ${LIBRARY}.testinputs )


ENDFOREACH( LIBRARY )

