
#
# This is a command file.
#
# To make a new test, all you have to do is:
#   1.  Make a new directory under tests/
#   2.  Put a file like this (named "command") into that directory.
#
# The contents of this file will be passed to the shell (Bash or SSH),
# so any legal shell commands can go in this file.
# Or comments like this one, for that matter.
#
# Variable substiution is done using Python's printf format,
# meaning you need a percent sign, the variable name in parentheses,
# and the letter 's' (for 'string').
#
# Available variables include:
#   workdir     the directory where test input files have been copied,
#               and where test output files should end up.
#   minidir     the base directory where Mini lives
#   database    where the Mini database lives
#   bin         where the Mini binaries live
#   binext      the extension on binary files, like ".linuxgccrelease"
#
# The most important thing is that the test execute in the right directory.
# This is especially true when we're using SSH to execute on other hosts.
# All command files should start with this line:
#

cd /Users/vmulligan/rosetta_devcopy/Rosetta/main/tests/integration/new/constraints_metric

[ -x /Users/vmulligan/rosetta_devcopy/Rosetta/main/source/bin/rosetta_scripts.default.macosclangdebug ] || exit 1

/Users/vmulligan/rosetta_devcopy/Rosetta/main/source/bin/rosetta_scripts.default.macosclangdebug  -database /Users/vmulligan/rosetta_devcopy/Rosetta/main/database -testing:INTEGRATION_TEST -info ConstraintsMetric 2>&1 \
    | egrep -vf ../../ignore_list \
    >info.log 

test "${PIPESTATUS[0]}" != '0' && exit 1 || true  # Check if the first executable in pipe line return error and exit with error code if so

/Users/vmulligan/rosetta_devcopy/Rosetta/main/source/bin/rosetta_scripts.default.macosclangdebug  -database /Users/vmulligan/rosetta_devcopy/Rosetta/main/database -parser:script_vars whichmetric=metric1 -testing:INTEGRATION_TEST -in:file:s xml/1l2y.pdb -in:file:fullatom -parser:protocol xml/test.xml 2>&1 \
    | egrep -vf ../../ignore_list \
    >out1.log 

test "${PIPESTATUS[0]}" != '0' && exit 1 || true  # Check if the first executable in pipe line return error and exit with error code if so

/Users/vmulligan/rosetta_devcopy/Rosetta/main/source/bin/rosetta_scripts.default.macosclangdebug  -database /Users/vmulligan/rosetta_devcopy/Rosetta/main/database -nstruct 2 -parser:script_vars whichmetric=metric2 -testing:INTEGRATION_TEST -in:file:s xml/1l2y.pdb -in:file:fullatom -parser:protocol xml/test.xml 2>&1 \
    | egrep -vf ../../ignore_list \
    >out2.log 

test "${PIPESTATUS[0]}" != '0' && exit 1 || true  # Check if the first executable in pipe line return error and exit with error code if so

test `grep CoordinateConstraint 1l2y_0001.pdb | wc -l` != '20' && exit 1 || true  # Check that the output has 20 full entries.
test `grep CoordinateConstraint 1l2y_0002.pdb | wc -l` != '10' && exit 1 || true  # Check that the partial output has 10 entries.