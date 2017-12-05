# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from __future__ import print_function
import timeit

import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.kinematics import FoldTree

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

pose1 = pyrosetta.pose_from_sequence('ACDEFGHI')

# Test __len__
assert(len(Pose().residues) == 0)
assert(len(pose1.residues) == 8)
assert(len(pose1) == 8) # Deprecated

# Test __iter__
for residue in Pose().residues:
    pass

for residue in pose1.residues:
    pass

for residue in pose1: # Deprecated
    pass


# Test __getitem__
# assert(pose1.residues[0] == ValueError)
# assert(pose1.residues[0:] == ValueError)
assert(pose1.residues[1].annotated_name() == 'A[ALA:NtermProteinFull]')
assert(pose1.residues[6].annotated_name() == 'G')
assert(pose1.residues[8].annotated_name() == 'I[ILE:CtermProteinFull]')
assert(pose1.residues[-1].annotated_name() == 'I[ILE:CtermProteinFull]')
assert(pose1.residues[-3].annotated_name() == 'G')
assert(pose1.residues[-8].annotated_name() == 'A[ALA:NtermProteinFull]')
assert(''.join([res.annotated_name() for res in Pose().residues[:]]) == '')
assert(''.join([res.annotated_name() for res in pose1.residues[:]]) == 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1.residues[1:9]]) == 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1.residues[:-3]]) == 'A[ALA:NtermProteinFull]CDEF')
assert(''.join([res.annotated_name() for res in pose1.residues[3:]]) == 'DEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1.residues[-6:8]]) == 'DEFGH')
assert(''.join([res.annotated_name() for res in pose1.residues[-6:8:2]]) == 'DFH')
assert(''.join([res.annotated_name() for res in pose1.residues[-6:8:3]]) == 'DG')
assert(''.join([res.annotated_name() for res in pose1[1:9]]) == 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]') # Deprecated
assert(''.join([res.annotated_name() for res in pose1[-6:8]]) == 'DEFGH') # Deprecated

# Test __iadd__
gly_residue = pose1.residues[6]
pose1.residues += gly_residue
assert(''.join([res.annotated_name() for res in pose1.residues]) == 'A[ALA:NtermProteinFull]CDEFGHIG[GLY:CtermProteinFull]')

pose2 = Pose()
pose2.residues += pose1.residues[1]
for _ in range(3):
    pose2.residues += gly_residue
pose2.residues += pose1.residues[-1]
pose2.residues += gly_residue
assert(''.join([res.annotated_name() for res in pose2.residues]) == 'A[ALA:NtermProteinFull]GGGGG[GLY:CtermProteinFull]')


# Test __imul__
pose3 = Pose()
pose3.residues *= pose1.residues[5]
assert(''.join([res.annotated_name() for res in pose3.residues]) == 'F')
pose3.residues *= pose1.residues[5]
assert(''.join([res.annotated_name() for res in pose3.residues]) == 'FF')
pose3.residues *= pose1.residues[5]
assert(''.join([res.annotated_name() for res in pose3.residues]) == 'FFF')

ft = FoldTree()
ft.add_edge(1, 2, 1)
ft.add_edge(2, 3, 2)
assert(pose3.fold_tree().__str__() == ft.__str__())

# Time Pose.residues
#print('Time to allocate 100000 residue accessors:',
#      timeit.timeit('pose_residues = pose.residues', setup='import pyrosetta; pose = pyrosetta.pose_from_sequence("ACDEFGHI")', number=100000))
