# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from __future__ import print_function
import timeit
from collections import defaultdict

import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.conformation import Residue
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector

pyrosetta.init(extra_options = "-constant_seed -out:mute all")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
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
# print('Time to allocate 100000 residue accessors:',
#        timeit.timeit('pose_residues = pose.residues', setup='import pyrosetta; pose = pyrosetta.pose_from_sequence("ACDEFGHI")', number=100000))

# Test Pose.residue_pairs
pose = pyrosetta.pose_from_file('../test/data/test_in_short.pdb')
for residue_pair in pose.residue_pairs():
    assert len(residue_pair) == 2
    assert all([isinstance(residue, Residue) for residue in residue_pair])

paired_residues = list()
for residue_pair in pose.residue_pairs(primary_residue_selector=ResidueIndexSelector(10)):
    assert(residue_pair[0].name1() == 'D')
    paired_residues.append(residue_pair[1].name1())
assert(len(paired_residues) == 19)
assert(''.join(paired_residues) == 'DAITIHSILWIEDNLESPL')

paired_residues = list()
for residue_pair in pose.residue_pairs(primary_residue_selector=ResidueIndexSelector(10), sequence_distance_minimum=5):
    assert(residue_pair[0].name1() == 'D')
    paired_residues.append(residue_pair[1].name1())
assert(len(paired_residues) == 9)
assert(''.join(paired_residues) == 'DAITLESPL')

paired_residues = list()
for residue_pair in pose.residue_pairs(secondary_residue_selector=ResidueIndexSelector(10), sequence_distance_minimum=1):
    assert(residue_pair[1].name1() == 'D')
    paired_residues.append(residue_pair[0].name1())
assert(len(paired_residues) == 17)
assert(''.join(paired_residues) == 'DAITIHSIIEDNLESPL')

paired_residues = list()
for residue_pair in pose.residue_pairs(primary_residue_selector=ResidueIndexSelector(10), neighborhood_distance_maximum=8.5):
    assert(residue_pair[0].name1() == 'D')
    paired_residues.append(residue_pair[1].name1())
assert(len(paired_residues) == 8)
assert(''.join(paired_residues) == 'HSILWIED')

paired_residues = list()
for residue_pair in pose.residue_pairs(primary_residue_selector=ResidueIndexSelector(10), neighborhood_distance_maximum=12.5):
    assert(residue_pair[0].name1() == 'D')
    paired_residues.append(residue_pair[1].name1())
assert(len(paired_residues) == 13)
assert(''.join(paired_residues) == 'ITIHSILWIEDNL')

paired_residues = list()
for residue_pair in pose.residue_pairs(primary_residue_selector=ResidueIndexSelector(10), sequence_distance_minimum=2, neighborhood_distance_maximum=8.5):
    assert(residue_pair[0].name1() == 'D')
    paired_residues.append(residue_pair[1].name1())
assert(len(paired_residues) == 4)
assert(''.join(paired_residues) == 'HSED')

paired_residues = defaultdict(list)
for residue_pair in pose.residue_pairs(primary_residue_selector=ResidueIndexSelector('5,10'), secondary_residue_selector=ResidueIndexSelector('5,6,7,8,9,10'), sequence_distance_minimum=2, neighborhood_distance_maximum=12.5):
    paired_residues[residue_pair[0].name1()].append(residue_pair[1].name1())

assert(len(paired_residues) == 2)
assert(paired_residues['D'] == ['I', 'H', 'S',])
assert(paired_residues['I'] == ['I', 'L', 'D'])

# Profile Pose.residue_pairs
# import cProfile
# pose = pyrosetta.pose_from_file("../test/data/test_in.pdb")
# def iter_residue_pairs():
#     for residue_pair in pose.residue_pairs(sequence_distance_minimum=2, neighborhood_distance_maximum=8.5):
#         pass
# cProfile.run('iter_residue_pairs()')
