#!/bin/sh
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
ctags          `find ../../src/core -name "*.cc" -or -name "*.hh"`
ctags --append `find ../../src/protocols -name "*.cc" -or -name "*.hh"`
ctags --append `find ../../src/devel -name "*.cc" -or -name "*.hh"`
ctags --append `find ../../src/utility -name "*.cc" -or -name "*.hh"`
ctags --append `find ../../src/numeric -name "*.cc" -or -name "*.hh"`
