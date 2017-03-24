# build test suite
These tests are mainly focused on verifying Rosetta build-ability in various modes.

-----
### debug
Compile Rosetta in debug mode using Scons. Test will pass if compilation is successful and fail otherwise.

-----
### release
Compile Rosetta release mode using Scons. Test will pass if compilation is successful and fail otherwise.

-----
### static
Compile Rosetta staticly in release mode using Scons. Test will pass if compilation is successful and fail otherwise.

-----
### ninja_debug
Compile Rosetta in debug mode using the CMake/Ninja build system. Test will pass if compilation is successful and fail otherwise.

-----
### ninja_release
Compile Rosetta in release mode using the CMake/Ninja build system. Test will pass if compilation is successful and fail otherwise.

-----
### ninja_graphics
Compile Rosetta's graphics build using the CMake/Ninja build system. Test will pass if compilation is successful and fail otherwise.

-----
### header
Check that all of Rosetta's headers each compile on their own. Test will pass if all headers compile individually and fail otherwise

-----
### levels
Check that there are no `#include`s of headers from incorrect library levels (e.g. no protocols headers in core).
Test will fail if there are any level-prohibited includes used.

-----
### cppcheck
Run the Cppcheck code quality analyzer on Rosetta. Test will fail if a code quality issue identified by Cppcheck is introduced.
Test will *continue* to fail for future tests until the issue is dealt with.
