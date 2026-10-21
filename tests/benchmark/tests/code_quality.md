# code quality test suite
These tests are focused on checking that Rosetta's codebase does not have any obvious issues.

-----

### serialization
Check if serialization routines are properly implemented for those classes which are serializable.

----
### clang_analysis
Run the Clang Static Analyzer on Rosetta and report any issues found.

----
### clang_tidy
Run Clang Tidy on Rosetta and report issues found.
(Does not modify code.)

-----
### cppcheck
Run the Cppcheck code quality analyzer on Rosetta. Test will fail if a code quality issue identified by Cppcheck is introduced.
Test will *continue* to fail for future tests until the issue is dealt with.

----
### beautification
Check if new code in branch if beautified according to Rosetta standards.

----
### beautify
Automatically beautify new (relative to the main branch) code in a branch and commit the results.

----
### submodule_regression
Checks to see if the current main commit has issues with the submodule versions.
In particular, it looks to see if it introduces a regression (removes submodule commits which are in rosetta's main branch)
or if the version trails behind the submodule's main branch.
