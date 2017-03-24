# unit test suite 
These tests are focused on running the unit test suite.

----
### (plain)
Run the unit test suite, comparing code behavior with a fixed reference behavior.
Fails if computed results do not match the (in-code) reference results.

-----
### valgrind
Run the Valgrind memory checking tool on the unit tests. 
Failure indicates there are memory issues with the failing tests.

-----
### valgrind_detailed
Run the Valgrind memory checking tool on the unit tests. 
Failure indicates there are memory issues with the failing tests.
Provides more diagnostic output than the plain valgrind run.

-----
### addsan
Run GCC's AddressSanitizer on the unit tests.
Failure indicates there are memory issues with the failing tests.

-----
### ubsan
Run GCC's UndefinedBehaviorSanitizer on the unit tests.
Failure indicates the failing tests invoke C++ undefined behavior and hence may have unstable/uncertain results.


