(c) Copyright Rosetta Commons Member Institutions.
(c) This file is part of the Rosetta software suite and is made available under license.
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
(c) For more information, see http://www.rosettacommons.org. Questions about this can be
(c) addressed to University of Washington CoMotion, email: license@uw.edu.

@author atom-moyer

Directory for applying persistent bindings to PyRosetta Rosetta objects. Activated upon `import pyrosetta`.

Important functions are found in `pyrosetta.bindings.utility`.

Example:

```
from pyrosetta.util import bind_method, bind_property
from pyrosetta.rosetta.core.pose import Pose

@bind_method(Pose)
def test_method(self):
    print('Now you can call pose.test_method()!')

@bind_property(Pose)
def test_property(self):
    return 'Now pose.test_property returns this'
```
