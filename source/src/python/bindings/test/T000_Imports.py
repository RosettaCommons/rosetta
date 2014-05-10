# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov

import rosetta
from rosetta import *

print Vector1([1, 2, 3, 4])
print Vector1([1, 2, 3, 4.5])
print Vector1(['a', 'b', 'c', 'd'])


print 'Testing Python bindings for Rosetta...'
print 'Init...'
rosetta.init()
print version()
