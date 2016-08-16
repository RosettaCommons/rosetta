# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import rosetta
import pyrosetta

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
