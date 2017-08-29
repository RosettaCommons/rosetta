# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# @author atom-moyer

Link to docs on dunder methods: https://docs.python.org/3/reference/datamodel.html#special-method-names

Dunder methods (Double UNDERscore methods) are built-in "magic" python methods
which define shortcut behavior specific to python syntax.  Dunder methods are
defined given preset method names which give them their name (ie __init__,
__add__, __len__... etc). For example, calling 2 + 2 is equivalent to
calling 2.__add__(2).  These shorthand notations allows for, simultaneously,
more concise and expressive code.

This directory is specifically set aside to bind (dunder) methods to the python
classes corresponding to the rosetta C++ classes such as Pose and Residue. These
are generally for dunder bindings which cannot be automatically generated like
__add__ or __getitem__ because each object requires implementing this
functionality uniquely and intelligently.

Practically, the system of binding the dunder methods is implemented via a
decorator function, bind_method. It is located in pyrosetta.util. This
decorator takes a previously defined class such as
pyrosetta.rosetta.core.pose.Pose and dynamically sets/patches a new attribute
to the class. In theory, any attribute could be patched onto a class, but this
is frowned upon. In PyRosetta, binding methods is limited to python centric
problems such as building the python API compatible with rosetta objects.

Example:
```
from pyrosetta.util import bind_method
from pyrosetta.rosetta.core.pose import Pose

@bind_method(Pose)
def __len__(self):
    return self.size()
```
