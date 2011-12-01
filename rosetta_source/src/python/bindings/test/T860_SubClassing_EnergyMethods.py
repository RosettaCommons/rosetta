#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov
## @brief  Demo for PyRosetta sub classing

import sys

import rosetta
import rosetta.core.scoring.methods

rosetta.init()

pose = rosetta.pose_from_pdb("test/data/test_in.pdb")


# rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy sub-classing -----------------------------------
@rosetta.EnergyMethod(version=2)  # version is optional here
class MyCI1B_Method(rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy):
    def __init__(self):
        rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy.__init__(self, self.creator() )

    def residue_energy(self, rsd, pose, emap):
        emap.get().set( self.scoreType, 2.0)

    # you can define this functions by hand. And if you don't @EnergyMethod will do that for you
    #def clone(self): return MyNewCI1B();
    #def version(self): return 141
    #def indicate_required_context_graphs(self, v): pass


sf_new = rosetta.core.scoring.ScoreFunction()
sf_new.set_weight(MyCI1B_Method.scoreType, 1)
print '---------------------------------------------'
print 'MyCI1B_Method Score:', sf_new.score(pose)


