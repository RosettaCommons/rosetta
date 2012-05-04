#!/usr/bin/env python
import itertools
import os
import subprocess
import shutil
import fnmatch
import re
import sys 
from optparse import OptionParser

mol_to_params = "~/rosetta/rosetta_source/src/python/apps/public/molfile_to_params.py" 

char_set = ['0','1','2','3','4','5',
			'6','7','8','9','A','B',
			'C','D','E','F','G','H',
			'I','J','K','L','M','N',
			'O','P','Q','R','S','T',
			'U','V','W','X','Y','Z']

def get_name_from_params(path,database):
	path=path.strip()
	#if path is commented, ignore it
	if len(path) == 0:
		return 0
	if path[0] == "#":
		return 0
	fullpath = database+path
	#print fullpath
	path = path.split("/")
	#is path a params file?
	type = path[-1].split(".").pop()
	if type == "params":
		if os.path.exists(fullpath):
			paramsfile = open(fullpath,'r')
			for line in paramsfile:
				line = line.split()
				if len(line) >0:
					if line[0] == "IO_STRING":
						return line[1]
		else:
			return 0
	else:
		return 0

def make_params(mol_path,ligand_name,base_name):
	logfile = "params/"+base_name+"/log.txt"
	log = open(logfile,'a')
	params_cmd = mol_to_params+" -n " +ligand_name + " " + mol_path
	params_child = subprocess.call(params_cmd,shell=True,stdout=log,stderr=log)
	log.close()
	if params_child is not 0:
		return params_child
	
	pdb_path = "params/"+base_name+"/"+ligand_name+"_conformers.pdb"
	pdb_file = open(pdb_path,'w')
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file,ligand_name+"_*.pdb"):
			conformer =open(file,'r')
			pdb_file.writelines(conformer.readlines())
			conformer.close()
			os.remove(file)
			
	pdb_file.close()
	
	if os.path.exists(ligand_name+".params"):
		paramsfile = open(ligand_name+".params",'a')
	else:
		return 1
	paramsfile.write("PDB_ROTAMERS "+ ligand_name+"_conformers.pdb\n")
	paramsfile.close()
	shutil.move(ligand_name+".params", "params/"+base_name+"/"+ligand_name+".params")
	return params_child
	
			

def get_disallowed_ligands(database):
	residue_type_set_path = database+"chemical/residue_type_sets/"
	residue_set_list = os.listdir(residue_type_set_path)
	disallowed_ligands = set()
	for residue_set in residue_set_list:
		if residue_set[0] != ".":
			current_residue_set_path = residue_type_set_path+residue_set+"/"
			residue_types = open(current_residue_set_path+"residue_types.txt",'r')
			for line in residue_types:
				name = get_name_from_params(line, current_residue_set_path)
				#print name
				if name != 0:
					disallowed_ligands.add(name)
	return disallowed_ligands


usage = "%prog -d path/to/rosetta/database --script_path=/path/to/molfile_to_params.py molfile_list.txt"
parser=OptionParser(usage)
parser.add_option("-d", dest ="database",help="path to minirosetta database",default=1)
parser.add_option("--script_path",dest="script",help="location of the molfile_to_params script",default = "")
(options, args) = parser.parse_args()

if not os.path.exists("params"):
	os.mkdir("params")

if options.script != "":
	mol_to_params = options.script


if not os.path.exists(mol_to_params):
	parser.error("ERROR: make sure to specify the path to the molfile_to_params script with the option --params")

if len(args) != 1:
	parser.error("ERROR: you must specify the path to a file containing a list of molfile paths")	

molfile_list_path = args[0]
ligand_names = itertools.product(char_set,repeat=3)

disallowed_ligands = get_disallowed_ligands(options.database+"/")


molfile_list = open(molfile_list_path,'r')
for molfile in molfile_list:
	molfile = molfile.strip()
	while True:
		ligand_name = ligand_names.next()
		ligand_name = "".join(ligand_name)
		if ligand_name not in disallowed_ligands:
			break
	mol_base_name = molfile.split("/").pop().split(".")[0]
	
	if os.path.exists("params/"+mol_base_name+"/"+ligand_name+".params"):
		print "ligand " + mol_base_name + " already processed, continuing"
		continue
	if not os.path.exists("params/"+mol_base_name):
		os.mkdir("params/"+mol_base_name)
	print ligand_name,mol_base_name
	params_status = make_params(molfile,ligand_name,mol_base_name)
	if params_status is not 0:
		print "WARNING: failed to generate params for " +mol_base_name
	
