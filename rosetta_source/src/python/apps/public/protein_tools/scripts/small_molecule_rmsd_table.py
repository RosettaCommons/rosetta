#!/usr/bin/env python2.5

import Bio.PDB
from rosettautil.protein import pdbStat
#from rosettautil.rosetta import rosettaScore
from rosettautil.util import fileutil
from optparse import OptionParser
import math
from numpy import linalg


def calculate_centroid(residue_atom_coords):
	atom_total=[0.0,0.0,0.0]
	number_atoms = 0
	for i in range(len(residue_atom_coords)):
		atom_total += residue_atom_coords[i]
#		print atom_total
		number_atoms += 1	
#	centroid = [x_total / len(residue_atom_coords),z_total / len(residue_atom_coords),z_total / len(residue_atom_coords)]
	centroid = [0.0,0.0,0.0]
	centroid = atom_total/number_atoms
#	print centroid
	return centroid 


def calculate_rms(model_1_atoms,model_2_atoms):
	assert(len(model_1_atoms)== len(model_2_atoms))
	d_2_sum = 0.0
	for i in range(len(model_1_atoms)):
		distance_2 = 0
		for j in [0,1,2]:
			distance_2+=(model_1_atoms[i][j]-model_2_atoms[i][j])**2
		d_2_sum += distance_2
	return math.sqrt(d_2_sum / len(model_1_atoms))    

#first we need to extract the atom coordinates of the option-defined res_name for model1 and model2

#then we can calculate the rmsd of the extracted atom coordinates
#	return pdbStat.atom_rms(native_atoms,decoy_atoms)

usage = "%prog [options] residue_name list of decoys"
parser=OptionParser(usage)
parser.add_option("--res_name",dest="residue_name",help="residue name of molecule you want to calculate rmsd of")
parser.add_option("--centroid",dest="centroid",help="calculate rmsd of molecule centroid",default=False,action="store_true")
parser.add_option("--full_atom",dest="full_atom",help="calculate rmsd of whole molecule")
parser.add_option("--include_hydrogen",dest="include_hydrogen",help="include hydrogen in rmsd calculation",default=False,action="store_true")
(options,args)= parser.parse_args()

if len(args) <2:
        parser.error("specify at least 2 protein to compare")

print "calculating RMSD of ", options.residue_name, " between ", len(args), " models"
if options.centroid:
	print "calculating centroid RMSD"
else:
	print "calculating all-atom RMSD"


#extract the coordinates of the files
list_model_coord=[]

number_list=[]

for model_file in args:
	model  = pdbStat.load_pdb(model_file)
#	print model_file
	atoms = model.get_atoms()
	residue_atom_coords = [] 
	number = model_file.split('_')[3].split('.')[0]
	number_list.append(number)
	for atom in atoms:
		res = atom.get_parent()
#		print res.get_resname()
		if res.get_resname() == options.residue_name:
#			print atom.get_name()
			if atom.get_name()[0] == "H" and not options.include_hydrogen:
				continue
#			print atom.get_name()[0]
			residue_atom_coords.append(atom.get_coord())

	list_model_coord.append(residue_atom_coords)
	print number,
print
#print coords_list[1]

x=len(args)
for i in range(0,x-1):
	print number_list[i],'0',
	if options.centroid:
		model_1_atoms = calculate_centroid(list_model_coord[i])
#		print model_1_atoms
	else:
		model_1_atoms = list_model_coord[i]
	for j in range(1,x):
		if j<=i:
			print '0',
		else:
			rms = 0
			if options.centroid:
				model_2_atoms = calculate_centroid(list_model_coord[j])
#				print "model_2_coords: ",model_2_atoms
				rms = linalg.norm(model_1_atoms-model_2_atoms)
			else:
				model_2_atoms = list_model_coord[j]
				rms = calculate_rms(model_1_atoms,model_2_atoms)
			print str(rms),
	print 	


