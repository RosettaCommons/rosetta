SET( LINK_TYPE SHARED )
### Any external libraries to link statically
SET( LINK_EXTERNAL_LIBS z dl stdc++ )

# Enable proper rpath handling on macos. This is required to set @rpath
# correctly in library names to enable relocatible linking.
set(CMAKE_MACOSX_RPATH 1)

# Use the build directory rpath during build linking, so that executables can
# be run from the build dir.
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# Set new rpaths during install. Add the target install prefix, linux relative
# paths, and macos relative paths to the install rpath.
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;\$ORIGIN/../lib;@executable_path/../lib")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH.
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
