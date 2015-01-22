DATABASE VALIDATION INTEGRATION TEST
====================================

This is an integration test to validate the database. Unlike most other integration tests,
this gives an absolute pass/fail for each test, rather than a failure due to differences.

Currently it only tests the validity of the Dunbrack library binaries.

Interpreting failures
---------------------

First check to see if there is an issue with your reference run, only then look at the
results in the new directory.

Error with "ref" directory
--------------------------

If you have an erroring test for you reference directory for a clean checkout
(Easy to see by greping "DOES NOT MATCH" for the log files in your ref directory),
you have corruption of your Dunbrack binary files in your Rosetta database.
Delete the *.bin files under $ROSETTA_DATABASE/rotamers/ and then re-run 
*WITH A CLEAN CHECKOUT* to force regeneration of the binary files.

Which binary files need to be deleted and regenerated should be listed in the log 
files that are showing the "DOES NOT MATCH" failures.

*DO NOT DELETE THE BINARY FILES UNLESS YOU'RE GETTING FAILURES WITH A *CLEAN* CHECKOUT.*
- Doing can mask potentially serious changes to the Rosetta codebase.

Error with "new" directory
--------------------------

If your ref directory is clean, but your new/ directory is showing "DOES NOT MATCH"
failures, this means you changed the binary reading/writing code. If intentional, you 
will need to update the binary version number in the current_binary_format_version_id() 
functions of src/core/pack/dunbrack/RotamerLibrary.cc. You may also need to change the 
operator==() functions of the SingleResidueDunbrackLibrary object and derived/related 
classes.
