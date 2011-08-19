#!/bin/sh
# This script is a hack to get Xcode to index the header files.
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

echo "// This file is a hack to get versions of Xcode prior to 3.1 to index headers" > headers.cc
find ../src ../test -name "*.hh" | sed -E 's/^(.+)$/#include "\1"/' >> headers.cc
