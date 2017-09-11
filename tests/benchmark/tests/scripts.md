# Scripts test suite
These tests are focused on running tests for Rosetta/PyRosetta scripts from rosetta_scripts and pyrosetta_scripts repositories.

-----
### rosetta
Run tests from rosetta_scripts repository. There are three tests for the rosetta_scripts_scripts repository; there are two main tests for the rosetta_scripts_scripts repository, and an auxiliary test.

### rosetta.validate
validate_all_scripts.py: tests a large number of XML files for whether they are valid according to Rosetta's internally-generated XML Schema. This script tests all of the .xml scripts that are listed in the rosetta_scripts_scripts/scripts_to_validate.py file. If a script requires any flags in order to be validated (most notably, the -parser:script_vars flag), then these flags may be specified by placing them in a flags file named X.flags to correspond X.xml.  (E.g., a script named ddmi.xml could have ddmi.flags). This .flags file should live in the same directory as the .xml file.

### rosetta.parse
parse_all_scripts.py: tests that the parse-my-tag methods of the objects defined in a large number of XML files can run without error (and also ensures that the XML files themselves are valid). This test is strictly more rigorous than the validate_all_scripts.py test; any file that would fail the validate_all_scripts test would also fail this test. This script tests all of the .xml scripts that are listed in the rosetta_scripts_scripts/sripts_to_parse.py file. If a script requires any flags in order to be parsed (e.g., the -parser:script_vars flag), then those flags may be specified by having placing them in a flags file named X.flags to correspond X.xml.  (E.g., a script named ddmi.xml could have ddmi.flags). This .flags file should live in the same directory as the .xml file.

### rosetta.verify
verify_all_scripts_accounted_for.py: tests that every .xml file that is in the repository is tested by either the first or the second test, or is in a list of .xmls that are meant to not be tested. This test ensures that people do not accidentally add files without deciding how they should be tested.
When you add an XML file to the repository, you should list it in one (and only one) of three files:
 -- rosetta_scripts_scripts/scripts_to_validate.py <--- for files to be tested by the validate_all_scripts.py test
 -- rosetta_scripts_scripts/scripts_to_parse.py <--- for files to be tested by the parse_all_scripts.py test
 -- rosetta_scripts_scripts/untested_scripts.py <--- for files that should be tested by neither script
There are four kinds of scripts that would be listed in the untested_scripts.py file:
 -- xmls that are xi:included into other scripts, and cannot independently be tested (because they are incomplete)
 -- scripts that are now out of date and that cannot be trivially updated with the rewriter,
 -- scripts that rely on code that has not yet been merged into master, and
 -- scripts that are symlinked to from the integration test directory as part of a full integration test and that don't need additional testing.

-----
### pyrosetta
Run tests from pyrosetta_scripts repository.
