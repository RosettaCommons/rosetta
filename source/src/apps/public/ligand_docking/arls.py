#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import sys

# By doing this before we even import the "real" script,
# we ensure reasonable error messages
# even with "syntax errors" like generator expressions.
if not hasattr(sys, "version_info") or sys.version_info < (2,3):
    print
    print "  *** Script requires Python 2.3 or higher! ***"
    print
    sys.exit(1)

if __name__ == "__main__":
    import arls_impl
    sys.exit(arls_impl.main(sys.argv[1:]))
