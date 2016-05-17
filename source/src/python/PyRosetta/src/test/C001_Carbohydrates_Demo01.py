#! /usr/bin/python
# This file simply contains a list of commands copied from a carbohydrates
# tutorial/demo that I am writing.  ~Labonte

from rosetta import *

init('-constant_seed -include_sugars')
import os; os.chdir('.test.output')

pm = PyMOL_Mover()


# Linear Oligosaccharides & IUPAC Sequences

from rosetta.core.pose import pose_from_saccharide_sequence

glucose = pose_from_saccharide_sequence('alpha-D-Glcp')
galactose = pose_from_saccharide_sequence('Galp')
mannose = pose_from_saccharide_sequence('->3)-a-D-Manp')
maltotriose = pose_from_saccharide_sequence('a-D-Glcp-' * 3)
isomaltose = pose_from_saccharide_sequence('->6)-Glcp-' * 2)
lactose = pose_from_saccharide_sequence('b-D-Galp-(1->4)-a-D-Glcp')

pm.apply(isomaltose)
pm.apply(glucose)
pm.apply(galactose)

print maltotriose
print isomaltose
print lactose

print maltotriose.chain_sequence(1)
print isomaltose.chain_sequence(1)
print lactose.chain_sequence(1)

for res in lactose: print res.seqpos(), res.name()
for res in maltotriose: print res.seqpos(), res.name()

print glucose.residue(1)
print galactose.residue(1)
print mannose.residue(1)


# Torsion Angles, PDB File HETNAM Records, & RingConformers

galactose.chi(1, 1)
galactose.chi(2, 1)
galactose.chi(3, 1)
galactose.chi(4, 1)
galactose.chi(5, 1)
galactose.chi(6, 1)

#observer = PyMOL_Observer(galactose, True)

galactose.set_chi(1, 1, 180)
galactose.set_chi(2, 1, 60)
galactose.set_chi(3, 1, 60)
galactose.set_chi(4, 1, 0)
galactose.set_chi(5, 1, 60)
galactose.set_chi(6, 1, -60)

maltotriose = pose_from_file('../test/data/carbohydrates/maltotriose.pdb')
isomaltose = pose_from_file('../test/data/carbohydrates/isomaltose.pdb')

pm.apply(maltotriose)

maltotriose.phi(1)
maltotriose.psi(1)
maltotriose.phi(1)
maltotriose.phi(2)
maltotriose.psi(2)
maltotriose.omega(2)
maltotriose.phi(3)
maltotriose.psi(3)

#observer = PyMOL_Observer(maltotriose, True)

for i in (2, 3):
    maltotriose.set_phi(i, 180)
    maltotriose.set_psi(i, 180)

pm.apply(isomaltose)

isomaltose.phi(2)
isomaltose.psi(2)
isomaltose.omega(2)

#observer = PyMOL_Observer(isomaltose, True)

isomaltose.set_phi(2, 180)
isomaltose.set_psi(2, 180)
isomaltose.set_omega(2, 180)

pm.apply(glucose)

Glc1 = glucose.residue(1)

for i in range(1, 6): print Glc1.nu(i)

print Glc1.ring_conformer(1)

ring_set = Glc1.type().ring_conformer_set(1)
conformer = ring_set.get_ideal_conformer_by_name('1C4')
glucose.set_ring_conformation(1, 1, conformer)
pm.apply(glucose)


# Modified Sugars, Branched Oligosaccharides, & PDB File LINK Records

LacNAc = pose_from_saccharide_sequence('b-D-Galp-(1->4)-a-D-GlcpNAc')
pm.apply(LacNAc)

Lex = pose_from_saccharide_sequence('b-D-Galp-(1->4)-[a-L-Fucp-(1->3)]-D-GlcpNAc')
pm.apply(Lex)

Lex = pose_from_file('../test/data/carbohydrates/Lex.pdb')
pm.apply(Lex)

print Lex
for i in range(2): print Lex.chain_sequence(i + 1)

for res in Lex: print res

Lex.phi(2)
Lex.psi(2)
Lex.phi(3)
Lex.psi(3)

Lex.chi(3, 1)
Lex.chi(4, 1)

Lex.chi(2, 1)
Lex.set_chi(2, 1, 180)
pm.apply(Lex)
Lex.chi(7, 1)
Lex.set_chi(7, 1, 0)
pm.apply(Lex)


# N- and O-Linked Glycans

N_linked = pose_from_file('../test/data/carbohydrates/N-linked_14-mer_glycan.pdb')
pm.apply(N_linked)
pm.send_ss(N_linked)
print N_linked
for i in range(4): print N_linked.chain_sequence(i + 1)

O_linked = pose_from_file('../test/data/carbohydrates/O_glycan.pdb')
pm.apply(O_linked)
pm.send_ss(O_linked)
print O_linked
for i in range(2): print O_linked.chain_sequence(i + 1)

N_linked.set_phi(N_linked.pdb_info().pdb2pose("B", 1), 180)
pm.apply(N_linked)
N_linked.set_psi(N_linked.pdb_info().pdb2pose("B", 1), 0)
pm.apply(N_linked)
N_linked.set_omega(N_linked.pdb_info().pdb2pose("B", 1), 90)
pm.apply(N_linked)

peptide = pose_from_sequence('ASA')
pm.apply(peptide)
from rosetta.core.pose.carbohydrates import glycosylate_pose, glycosylate_pose_by_file
glycosylate_pose(peptide, 2, 'Glcp')
pm.apply(peptide)

peptide = pose_from_sequence('ASA')
glycosylate_pose_by_file(peptide, 2, '../database/chemical/carbohydrates/common_glycans/core_5_O-glycan.iupac')
pm.apply(peptide)
