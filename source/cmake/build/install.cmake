# Disable implicit dependency on "all" target for install.
# Allows installation of limited set of executables by indicating
# them as the primary build target.
# Eg. ninja rosetta_scripts score relax && ninja install
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

# Setup target destination for external and rosetta libraries.
# cppdb and it's support libraries (INTERNAL_LIBRARIES) are defined
# in external.cmake and aren't added to EXTERNAL_LIBRARIES
install(
  TARGETS ${LIBRARIES} ${EXTERNAL_LIBRARIES} cppdb ${INTERNAL_LIBRARIES}
  LIBRARY DESTINATION lib
)

# Copy all database files into the install directory.
install(
  DIRECTORY ../../../database/
  DESTINATION database
)

# Use score to trigger a build of any needed database files, ie binary dunbrack files.
install(CODE "MESSAGE(\"Executing score to setup rosetta database.\")")
install(CODE "MESSAGE(\"${CMAKE_INSTALL_PREFIX}/bin/score -out:file:scorefile /dev/null -in:file:s ${CMAKE_BINARY_DIR}/../../test/core/pack/dunbrack/1UBQ_repack.pdb\")")
install(CODE "execute_process(COMMAND ${CMAKE_INSTALL_PREFIX}/bin/score -out:file:scorefile /dev/null -in:file:s ${CMAKE_BINARY_DIR}/../../test/core/pack/dunbrack/1UBQ_repack.pdb)")
