# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta
import pyrosetta.rosetta as rosetta

pyrosetta.init()
print( pyrosetta.version() )

# Testing various ReturnValuePolicy

print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_PoseOP() )
print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_PoseCOP() )
# print rosetta.protocols.toolbox.PyReturnValuePolicyTest_PoseAP()
# print rosetta.protocols.toolbox.PyReturnValuePolicyTest_PoseCAP()


print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_ScoreFunctionOP() )
print( 'rosetta.protocols.toolbox.PyReturnValuePolicyTest_ScoreFunctionOP() PASSED' )

print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_ScoreFunctionCOP() )
print( 'rosetta.protocols.toolbox.PyReturnValuePolicyTest_ScoreFunctionCOP() PASSED' )

print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_ScoreFunctionCOP2() )
print( 'rosetta.protocols.toolbox.PyReturnValuePolicyTest_ScoreFunctionCOP2() PASSED' )


print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_SF_ReplicaOP() )
print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_SF_ReplicaCOP() )
# print rosetta.protocols.toolbox.PyReturnValuePolicyTest_SF_ReplicaAP()
# print rosetta.protocols.toolbox.PyReturnValuePolicyTest_SF_ReplicaCAP()


print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_DummyClassOP() )
print( rosetta.protocols.toolbox.PyReturnValuePolicyTest_DummyClassCOP() )
# print rosetta.protocols.toolbox.PyReturnValuePolicyTest_DummyClassAP()
# print rosetta.protocols.toolbox.PyReturnValuePolicyTest_DummyClassCAP()


from rosetta.protocols.toolbox.py_inheritance_test import *

base   = Base()

o_1A = O_1A()
o_1B = O_1B()

o_2A = O_2A()
o_2B = O_2B()

o_2A2B = O_2A2B()


def py_take_Base(o):
    take_Base_reference(o)
    take_Base_pointer(o)
    take_Base_SP(o)


def py_take_O_1A(o):
    take_O_1A_reference(o)
    take_O_1A_pointer(o)
    take_O_1A_SP(o)


def py_take_O_2A(o):
    take_O_2A_reference(o)
    take_O_2A_pointer(o)
    take_O_2A_SP(o)


for o in base, o_1A, o_1B, o_2A, o_2B, o_2A2B: py_take_Base(o)

for o in       o_1A,       o_2A,       o_2A2B: py_take_O_1A(o)

for o in                   o_2A,       o_2A2B: py_take_O_2A(o)
