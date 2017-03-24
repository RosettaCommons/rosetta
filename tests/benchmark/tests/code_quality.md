# code quality test suite
These tests are focused on checking that Rosetta's codebase does not have obvious issues.

-----

### serialization
Check if serialization routines are properly implemented for those classes which are serializable.

----
### clang_analysis
Run the Clang Static Analyzer on Rosetta and report any issues found.

----
### beautification
Check if new code in branch if beautified according to Rosetta standards.

----
### beautify
Automatically beautify new (relative to master) code in a branch and commit the results.
