import sys
import PSSM
from Bio.PDB import * 
from rosettautil.util import fileutil
import math

def load_pdb(path):
	"""return a biopython structure object given a pdb file path"""
	parser = PDBParser(PERMISSIVE=1)
	pdb_file = fileutil.universal_open(path,'rU')
	structure = parser.get_structure(path[0:4],pdb_file)
	pdb_file.close()
	return structure