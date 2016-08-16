#!/bin/bash
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#UTracers is a system for doing file output comparisons in unit
#tests. Note: Consider making an integration test instead.

#This script replaces the generated UTracers with the reference ones
#stored in source/test Note: if you use non-standard build with extras
#you'll probably need to modify the sed statement below.

for i in $(find  build -name "*._tmp_"); do cp "$i" $( echo ${i%%._tmp_} | sed 's/.*\/default\//test\//' ); done

# scons--in its infinite wisdom--does not recognized changes to
# UTracers. So, wipe the test directory to force a rebuild.
rm -rf build/test
