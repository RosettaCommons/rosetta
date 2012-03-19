import sys
import PSSM
from Bio.PDB import * 
from Bio import pairwise2
from rosettautil.util import fileutil
import math
import warnings


def sequence_recovery(native_struct,designed_struct):
	"""calculate percent sequence recovery between a native and designed struct"""
	native_residues = native_struct.get_residues()
	designed_residues = designed_struct.get_residues()
	total = 0.0;
	recovered = 0.0;
	for native,designed in zip(native_residues,designed_residues):
		if native.get_resname() == designed.get_resname():
			recovered += 1
		total += 1
	#print recovered, total
	return recovered/total

def sequence_recovery_range(native_struct,designed_struct,min,max):
	"""calculate percent sequence recovery between a native and designed
	struct for all residues between a minimum and maximum b factor """
	native_residues = native_struct.get_residues()
	designed_residues = designed_struct.get_residues()
	total = 0.0
	recovered = 0.0
	for native,designed in zip(native_residues,designed_residues):
		score = native.get_list()[1].get_bfactor()
		if score >= min and score <= max:
			if native.get_resname() == designed.get_resname():
				recovered += 1
			total += 1
	return recovered/total

def sequence_recovery_group(native_struct,designed_struct,min,max):
	"""calculate sequence recovery by group (non-polar, polar,
	aromatic, and charged) for all residues between a minimum and
	maximum b factor"""
	non_polar = ["GLY","ALA","VAL","LEU","MET","ILE"]
	polar = ["SER","THR","CYS","PRO","ASN","GLN"]
	aromatic = ["PHE","TYR","TRP"]
	charged = ["LYS","ARG","HIS","ASP","GLU"]

	groups = {"non_polar":non_polar,"polar":polar,"aromatic":aromatic,"charged":charged}
	

	native_residues = native_struct.get_residues()
	designed_residues = designed_struct.get_residues()
	total = {"non_polar" :0.0, "polar":0.0, "aromatic":0.0, "charged":0.0}
	recovered = {"non_polar" :0.0, "polar":0.0, "aromatic":0.0, "charged":0.0}

	for native,designed in zip(native_residues,designed_residues):
		score = native.get_list()[1].get_bfactor()
		
		if score >= min and score <= max:
			for group in groups:
				if native.get_resname() in groups[group] and designed.get_resname() in groups[group]:
					recovered[group] += 1
				if designed.get_resname() in groups[group]:
					total[group] += 1
					break
	for group in groups:
		print group,recovered[group],total[group]
		#recovered[group]=recovered[group]/total[group]

	return recovered

	

def sequence_composition(struct):
	"""calculate sequence composition by residue"""
	struct_residues = struct.get_residues()
	composition = {}
	for residue in struct_residues:
		residue_name = residue.get_resname()
		try:
			composition[residue_name]+=1
		except KeyError:
			composition[residue_name] = 1 
	return composition
	
def sequence_composition_range(struct,min,max):
	"""calculate sequence composition for all residues within a
	minimum and maximum b factor"""
	struct_residues = struct.get_residues()
	composition = {}
	for residue in struct_residues:
		score = residue.get_list()[1].get_bfactor()
		if score >= min and score <= max:
			residue_name = residue.get_resname()
			try:
				composition[residue_name]+=1
			except KeyError:
				composition[residue_name] = 1 
	return composition

def pssm_recovery_map(struct,pssm_map):
	"""calculate the pssm recovery given a structure and a pssm map"""
	struct_residues = struct.get_residues()
	#pssm_recovery = 0.0;
	#struct_size = 0.0;
	recovery_map = {}
	for residue in struct_residues:
		residue_name = residue.get_resname()
		residue_num = residue.get_id()[1]
		status = pssm_map.conserved(residue_num,residue_name)
		if status:
			try:
				recovery_map[residue_name]+=1
			except KeyError:
				recovery_map[residue_name] = 1
	return recovery_map

def pssm_recovery_map_range(struct,pssm_map,min,max):
	"""calculate the pssm recovery within a range of b factors given a
	structure and a pssm map"""
	struct_residues = struct.get_residues()
	recovery_map = {}
	for residue in struct_residues:
		score = residue.get_list()[1].get_bfactor()
		if score >= min and score <= max:
			residue_name = residue.get_resname()
			residue_num = residue.get_id()[1]
			status = pssm_map.conserved(residue_num, residue_name)
			if status:
				try:
					recovery_map[residue_name]+= 1
				except KeyError:
					recovery_map[residue_name] = 1
	return recovery_map

def pssm_recovery(struct,pssm_map):
	"""return percent pssm recovery"""
	struct_residues = struct.get_residues()
	pssm_recovery = 0.0;
	struct_size = 0.0;
	for residue in struct_residues:
		residue_name = residue.get_resname()
		residue_num = residue.get_id()[1]
		status = pssm_map.conserved(residue_num,residue_name)
		if status:
			pssm_recovery += 1.0
		struct_size += 1.0
	return pssm_recovery/struct_size

def pssm_recovery_range(struct,pssm_map,min,max):
	"""return percent pssm recovery fro residues within a range of b factors"""
	pssm_recovery = 0.0;
	struct_size = 0.0;
	for residue in struct.get_residues():
		score= residue.get_list()[1].get_bfactor()
		#print score
		if score >= min and score <= max:
			residue_name = residue.get_resname()
			residue_num = residue.get_id()[1]
			status = pssm_map.conserved(residue_num,residue_name)
			if status:
				pssm_recovery += 1.0
			struct_size += 1.0
			
	return pssm_recovery/struct_size

def pssm_scores(struct,pssm_map):
	""" return raw pssm total for each residue"""
	struct_residues = struct.get_residues()
	pssm_scores = {}
	size = 0
	for residue in struct_residues:
		size += 1
		residue_name = residue.get_resname()
		residue_num = residue.get_id()[1]
		try:
			pssm_scores[residue_name] += pssm_map.get_score(residue_num,residue_name)
		except KeyError:
			pssm_scores[residue_name] = pssm_map.get_score(residue_num,residue_name)
	for key in pssm_scores:
		pssm_scores[key] = pssm_scores[key]/size
	return pssm_scores


def pssm_scores_range(struct,pssm_map,min,max):
	""" return raw pssm total for each residue within a range of pssm scores"""
	struct_residues = struct.get_residues()
	pssm_scores = {}
	size = 0
	for residue in struct_residues:
		score = residue.get_list()[1].get_bfactor()
		#print score
		if score >= min and score <= max:
			size += 1
			residue_name = residue.get_resname()
			residue_num = residue.get_id()[1]
			try:
				pssm_scores[residue_name] += pssm_map.get_score(residue_num,residue_name)
			except KeyError:
				pssm_scores[residue_name] = pssm_map.get_score(residue_num,residue_name)
	for key in pssm_scores:
		pssm_scores[key] = pssm_scores[key]/size
	return pssm_scores


def ca_rms_only(struct_a,struct_b,residue_list):
	"""calculate the CA rmsd of two structures using the residues in residue_list"""
	residues_a = struct_a.get_residues();
	residues_b = struct_b.get_residues();

	d_2_sum = 0.0
	resn = 0
	for (res_a, res_b) in zip(residues_a, residues_b):
		if res_a.get_id()[1] not in residue_list:
			continue
		CA_a = res_a['CA']
		CA_b = res_b['CA']

		distance_2 = (CA_a-CA_b)**2
		d_2_sum += distance_2
		resn += 1
		
	rmsd = math.sqrt(d_2_sum/resn)
	return rmsd


def atom_rms(atoms_a,atoms_b,residue_list=""):
	"""calculate the all atom rmsd given two lists of atoms and a list of residues"""
	d_2_sum = 0.0
	resn = 0
	for(atom_a,atom_b) in zip(atoms_a,atoms_b):
		parent_a = atom_a.get_parent()
		if parent_a.get_id()[1] in residue_list or residue_list == "":
			#print "calculating for",parent_a.get_id()[1]
			distance_2 = (atom_a-atom_b)**2
			d_2_sum += distance_2
			resn +=1
		else:
			continue
	rmsd = math.sqrt(d_2_sum/resn)
	return rmsd	

def copy_b_factor(native_pdb,designed_pdb):
	""" copy the b factors from one structure to another"""
	native_atoms = native_pdb.get_atoms()
	designed_atoms = designed_pdb.get_atoms()
	for native_atom, designed_atom in zip(native_atoms, designed_atoms):
		designed_atom.set_bfactor(native_atom.get_bfactor)
	return designed_pdb

def calculate_rms(native,decoy,ca_mode,residues,rms_residues,chain):
    native_atoms = []
    decoy_atoms = []
    residue_set = set()
    rms_residue_set = set()
    if residues is not "":
        residue_file = open(residues,'r')
        for residue in residue_file:
            residue=residue.strip()
            if residue is not '':
                residue_set.add(int(residue))
        residue_file.close()

    if rms_residues is not "":
        rms_residue_file = open(rms_residues,'r')
        for rms_residue in rms_residue_file:
            rms_residue = rms_residue.strip()
            if rms_residue is not '':
                rms_residue_set.add(int(rms_residue))
        rms_residue_file.close()
        
    for(native_chain, decoy_chain) in zip (native[0].get_list(),decoy[0].get_list()):
        if chain is not "":
           if native_chain.get_id() != chain:
               print "ignoring chain " + native_chain.get_id()
               continue
        for(native_residue, decoy_residue) in zip(native_chain.get_list(),decoy_chain.get_list()):
           # if rms_residues is not "":
           #     if native_residue.get_id()[1] not in rms_residue_set:
           #         continue
            if len(residue_set) > 0 and native_residue.id[1] in residue_set:
                if ca_mode:
                    try:
                        native_atoms.append(native_residue['CA'])
                        decoy_atoms.append(decoy_residue['CA'])
                    except KeyError:
                        print "WARNING: residue", str(native_residue.get_id()[1]), "has no CA atom in either the native or the decoy structure.  Either this residue is a hetatm or part of your backbone is missing. The residue will not be included in the RMSD calculations."
                else:
                    for (native_atom,decoy_atom) in zip(native_residue.get_list(), decoy_residue.get_list()):
                       # print native_atom
			native_atoms.append(native_atom)
                        decoy_atoms.append(decoy_atom)
            elif len(residue_set) is 0:
                if ca_mode:
                    try:
                        native_atoms.append(native_residue['CA'])
                        decoy_atoms.append(decoy_residue['CA'])
                    except KeyError:
                        print "WARNING: residue", str(native_residue.get_id()[1]), "has no CA atom in either the native or the decoy structure.  Either this residue is a hetatm or part of your backbone is missing. The residue will not be included in the RMSD calculations."
                else:
                    for (native_atom,decoy_atom) in zip(native_residue.get_list(), decoy_residue.get_list()):
                        #print native_atom
			native_atoms.append(native_atom)
                        decoy_atoms.append(decoy_atom)
    superpose = Superimposer()
    superpose.set_atoms(native_atoms,decoy_atoms)
    superpose.apply(decoy.get_atoms())
    if rms_residues is not "":
        return pdbStat.atom_rms(native_atoms,decoy_atoms,rms_residue_set)
    else:
        return superpose.rms

def find_gaps(pdb,sequence,chain):
    """return a sequence with gaps in the pdb represented by - symbols"""
    #we can't trust the seqres record, it might not even exist, so get the sequence by looping through all the residues
    pdb_sequence = ""
    for residue in pdb.get_residues():
        if residue.get_full_id()[2] != chain: #skip residues that aren't in the chain
            continue
        residue_name = Polypeptide.three_to_one(residue.get_resname())
        type(residue_name)
        pdb_sequence += residue_name
    #now we align the two sequences
    alignment = pairwise2.align.globalmx(pdb_sequence,sequence,2,-1)
    #print alignment
    return alignment[0][1]
    
        
