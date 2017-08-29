# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

pose1 = pyrosetta.pose_from_sequence('ACDEFGHI')

# Test __len__
assert(len(pose1) == 8)

# Test __iter__
for residue in pose1:
    pass

# Test __getitem__
assert(pose1[1].annotated_name() == 'A[ALA:NtermProteinFull]')
assert(pose1[6].annotated_name() == 'G')
assert(pose1[8].annotated_name() == 'I[ILE:CtermProteinFull]')
assert(pose1[-1].annotated_name() == 'I[ILE:CtermProteinFull]')
assert(pose1[-3].annotated_name() == 'G')
assert(pose1[-8].annotated_name() == 'A[ALA:NtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1[:]]) == 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1[0:9]]) == 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1[1:9]]) == 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1[:-3]]) == 'A[ALA:NtermProteinFull]CDEF')
assert(''.join([res.annotated_name() for res in pose1[3:]]) == 'DEFGHI[ILE:CtermProteinFull]')
assert(''.join([res.annotated_name() for res in pose1[-6:8]]) == 'DEFGH')
assert(''.join([res.annotated_name() for res in pose1[-6:8:2]]) == 'DFH')
assert(''.join([res.annotated_name() for res in pose1[-6:8:3]]) == 'DG')

# Test __iadd__
gly_residue = pose1[6]
pose1 += gly_residue
assert(''.join([res.annotated_name() for res in pose1]) == 'A[ALA:NtermProteinFull]CDEFGHIG[GLY:CtermProteinFull]')

pose2 = pyrosetta.Pose()
pose2 += pose1[1]
for _ in range(3):
    pose2 += gly_residue
pose2 += pose1[-1]
pose2 += gly_residue
assert(''.join([res.annotated_name() for res in pose2]) == 'A[ALA:NtermProteinFull]GGGGG[GLY:CtermProteinFull]')

# Test __imul__
pose3 = pyrosetta.Pose()
pose3 *= pose1[5]
assert(''.join([res.annotated_name() for res in pose3]) == 'F')
pose3 *= pose1[5]
assert(''.join([res.annotated_name() for res in pose3]) == 'FF')
