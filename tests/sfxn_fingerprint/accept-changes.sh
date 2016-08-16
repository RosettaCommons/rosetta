#!/bin/bash
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


# Simple shortcut to copy all test results from new/ into ref/,
# thereby updating the "expected" test results.

for d in $( ls -d new/*); do rm -r ref/$(basename $d); mv -v $d ref/; done
