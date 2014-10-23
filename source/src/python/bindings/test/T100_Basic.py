# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov

from rosetta import *

rosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


print 'General testing ----------------------------------------------'

print 'Creating Pose object...'
pose = Pose()

print 'Pose from PDB...'
pose = pose_from_pdb("../test/data/test_in.pdb")

# TODO: rename pose_from_pdb or make_pose_from_sequence to be parallel

print 'Building Pose from sequence...'
pose3 = Pose()
make_pose_from_sequence(pose3, "DSEEKFLRRIGRFGYGYGPYE",'centroid')
print pose3

pose4 = Pose()
make_pose_from_sequence(pose4, "ARNDCEQGHILKMFPSTWYV", 'fa_standard')


print 'Dump PDB...'
dump_pdb(pose, "T110_Basic._.pdb")

print 'accessing pose attributes'
print pose
# TODO: remove extra blank lines at end
print 'there are ', pose.total_residue(), 'residues in this pose object'
print 'phi of residue 5 is ', pose.phi(5)
print 'psi of residue 5 is ', pose.psi(5)

print 'set phi of residue 5 to -60'
pose.set_phi(1, -60)
print 'set psi of residue 5 to -50'
pose.set_psi(1, -50)

print 'accessing residue 5 from pose'
res5 = pose.residue(5)
print res5

print 'accessing atoms from residue 5'
at5N  = res5.atom('N')
at5CA = res5.atom("CA")
at5C  = res5.atom("C")

print at5N
# TODO: above should print atom type key not magic number
# 2/23/9: hm, not sure this is possible b/c atom does not know which AtomTypeSet to use!

print 'xyz of at5N:', at5N.xyz().x, at5N.xyz().y, at5N.xyz().z
print 'norm of xyz at5N:', at5N.xyz().norm

print res5.atoms()  # <-- Still missing

atomN = AtomID(1,5)
atomCA = AtomID(2,5)
atomC = AtomID(3,5)
print 'bond length of N-CA in residue 5 is '
print pose.conformation().bond_length(atomN, atomCA)
print 'bond angle of N-CA-C in residue 5 is '
print pose.conformation().bond_angle(atomN, atomCA, atomC)
print 'setting bond length of N-CA in residue 5 to 1.5A '
pose.conformation().set_bond_length(atomN, atomCA, 1.5)
print 'setting bond angle of N-CA-C in residue 5 to 90 '
pose.conformation().set_bond_angle(atomN, atomCA, atomC, 90)
# TODO: make the above work with atom objects instead of atomIDs


print 'pose was generated from this pdb file: ', pose.pdb_info().name()
print 'pose numbering for chain A, residue 5, is ', pose.pdb_info().pdb2pose('A',5)
print 'pdb chain letter and residue number for residue 5, is ', pose.pdb_info().pose2pdb(5)
# TODO: pdb_info.* does not tab-complete


# Creating residue example
chm = rosetta.core.chemical.ChemicalManager.get_instance()
rts = chm.residue_type_set('fa_standard')
ala = rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('ALA') )
print ala
