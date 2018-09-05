#! /usr/bin/python
# List of commands used in PyRosetts Workshop #2

from __future__ import print_function

from math import *

# Basic PyRosetta
from rosetta import *
from pyrosetta import *
from pyrosetta.toolbox import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


pose = pose_from_file("../test/data/workshops/1YY8.clean.pdb")

# Commenting out this for now to allow tests pass without net connection
#pose2 = pose_from_rcsb("1YY8")


print( pose )
print( pose.sequence() )
print( "Protein has", pose.total_residue(), "residues." )
print( pose.residue(500).name() )

print( pose.pdb_info().chain(500) )
print( pose.pdb_info().number(500) )

print( pose.pdb_info().pdb2pose("A", 100) )

print( pose.pdb_info().pose2pdb(25) )

pose.display_secstruct()

# Protein Geometry
print( pose.phi(5) )
print( pose.psi(5) )
print( pose.chi(1, 5) )
R5N = AtomID(1, 5)
R5CA = AtomID(2, 5)
R5C = AtomID(3, 5)
print( pose.conformation().bond_length(R5N, R5CA) )
print( pose.conformation().bond_length(R5CA, R5C) )
N_xyz = pose.residue(5).xyz("N")
CA_xyz = pose.residue(5).xyz("CA")
N_CA_vector = CA_xyz - N_xyz
print( N_CA_vector.norm )

print( pose.conformation().bond_angle(R5N, R5CA, R5C) )
pose.set_phi(5, -60)
pose.set_psi(5, -43)
pose.set_chi(1, 5, 180)

pose.conformation().set_bond_length(R5N, R5CA, 1.5)
pose.conformation().set_bond_angle(R5N, R5CA, R5C, 110./180.*3.14159)
