#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov
## @brief  Demo for PyRosetta sub classing

from __future__ import print_function

import sys

import rosetta, pyrosetta
#import rosetta.utility.py
#import rosetta.core.pose
#import rosetta.core.scoring
#import rosetta.core.scoring.methods

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


class PyValue(rosetta.utility.py.Value):
    def __init__(self, v='PyUnknown'):
        print( 'PyValue.__init__...' )
        rosetta.utility.py.Value.__init__(self, v)


# Dummy C++ class sub-classing
class A(rosetta.utility.py.Base):
    def __init__(self):
        print( 'A.__init__...' )
        rosetta.utility.py.Base.__init__(self)

    def foo_int(self, i):    print( "A.foo_int({})".format(i) )
    def foo_string(self, s): print( "A.foo_string({})".format(s) )

    def foo_value(self, v):     print( "A.foo_value({})".format(v.s()) );    v.s('A.foo_value!')
    def foo_value_sp(self, v):  print( "A.foo_value_sp({})".format(v.s()) ); v.s('A.foo_value_sp!')


a = A()
v = rosetta.utility.py.Value()
pv = PyValue()

v.s('wrong...')
rosetta.utility.py.foo_all(a, 1, v, "qwerty" )
print( 'New Value[1]:', v.s(), '\n' )
assert v.s() == 'A.foo_value!'

pv.s('wrong...')
a.foo_value(pv)
print( 'New Value[2]:', pv.s(), '\n' )
assert pv.s() == 'A.foo_value!'

v.s('wrong...')
a.foo_value_sp(v)
print( 'New Value[3]:', v.s(), '\n' )
assert v.s() == 'A.foo_value_sp!'

pv.s('wrong...')
rosetta.utility.py.foo_all(a, 1, pv, "qwerty" )
print( 'New Value[4]:', pv.s(), '\n' )
assert pv.s() == 'A.foo_value!'

pv.s('wrong...')
rosetta.utility.py.foo_all_sp(a, 1, pv, "qwerty" )
print( 'New Value[5]:', pv.s(), '\n' )
assert pv.s() == 'A.foo_value!'


# Dummy C++ class sub-classing
class PyOverloadTest(rosetta.utility.py.OverloadTest):
    def __init__(self):
        rosetta.utility.py.OverloadTest.__init__(self)

    def test_p(self, v): print('PyOverloadTest.test_p({})'.format(v))
    def pure_test_p(self, v): print('PyOverloadTest.pure_test_p({})'.format(v))

    def test_ref(self, v): print('PyOverloadTest.test_ref({})'.format(v))
    def pure_test_ref(self, v): print('PyOverloadTest.pure_test_ref({})'.format(v))

print('rosetta.utility.py.OverloadTest Test -----------------------')
ot = rosetta.utility.py.OverloadTest()
ot.self_test_virtual_p()
ot.self_test_virtual_ref()

print('PyOverloadTest Test -----------------------------------------')
pot = PyOverloadTest()
pot.self_test_virtual_p()
pot.self_pure_test_virtual_p()
pot.self_test_virtual_ref()
pot.self_pure_test_virtual_ref()

# Mover sub-classing -----------------------------------
class My_New_Mover(rosetta.protocols.moves.Mover):
    def __init__(self):
        print( 'My_New_Mover.__init__...' )
        rosetta.protocols.moves.Mover.__init__(self)

    def get_name(self): return 'My_New_Mover'

    def apply(self, p):
        print( 'My_New_Mover.apply:', type(p) )
        # if isinstance(p, rosetta.core.pose.PoseAP):
        #     print 'Got rosetta.core.pose.PoseAP!'
        #     p = p.get()

        print( 'This My_New_Mover apply...' )
        p.set_phi(1, p.phi(1)+1.54)


new_mover = My_New_Mover()

pose = rosetta.core.import_pose.pose_from_file("../test/data/test_in.pdb")
sf_new = rosetta.core.scoring.ScoreFunction()

#new_mover.apply(pose)

seq = rosetta.protocols.moves.SequenceMover()
seq.add_mover( new_mover )

minmover = rosetta.protocols.simple_moves.MinMover()

old_phi = pose.phi(1)
seq.apply(pose)
print( 'Old phi=%s, new phi=%s' % (old_phi, pose.phi(1)) )
if pose.phi(1) == old_phi : sys.exit(1)


# rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy sub-classing -----------------------------------

# Making energy creator class by hand, this is for demo purpose only, in real programm use class decorator instead (for details see T860_SubClassing_EnergyMethods.py and T870_SubClassing_EnergyMethods2.py)
_mem_ = []
class MyNewCI1B_Creator(rosetta.core.scoring.methods.EnergyMethodCreator):
    def __init__(self):
        rosetta.core.scoring.methods.EnergyMethodCreator.__init__(self)

    def create_energy_method(self, energy_method_options):
        e = MyNewCI1B()
        _mem_.append(e)
        return e

    def score_types_for_method(self):
        sts = rosetta.utility.vector1_core_scoring_ScoreType();  sts.append( rosetta.core.scoring.PyRosettaEnergy_last )
        return sts


class MyNewCI1B(rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy):
    def __init__(self):
        print( 'MyNewCI1B::__init__!' )
        rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy.__init__(self, MyNewCI1B_Creator() )

    def residue_energy(self, rsd, pose, emap):
        #print( 'residue_energy:', type(emap) )
        emap.set(rosetta.core.scoring.PyRosettaEnergy_last, 2.0)


    def clone(self): return MyNewCI1B();

    def version(self): return 1

    def indicate_required_context_graphs(self, v): pass


b = MyNewCI1B_Creator()
rosetta.core.scoring.methods.PyEnergyMethodRegistrator( b )


sf_new = rosetta.core.scoring.ScoreFunction()
sf_new.set_weight(rosetta.core.scoring.PyRosettaEnergy_last, 1)
print( '---------------------------------------------' )
print( 'Score:', sf_new.score(pose) )
assert sf_new.score(pose) == 232.0
