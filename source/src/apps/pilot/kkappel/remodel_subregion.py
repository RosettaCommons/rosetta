import os
import argparse
from get_sequence import get_sequences
from read_pdb import read_pdb
from make_tag import make_tag_with_dashes
import rna_helper

# a script to setup density modeling jobs for sub-regions

def get_coord_csts( ref_pdb, exclude_resi=[] ):

	#func = 'FADE -10.0 10.0 5.0 -100.0 100.0' # no penalty if positions deviate up to 5A from initial coords
	func = 'FADE -15.0 15.0 5.0 -100.0 100.0' # no penalty if positions deviate up to 5A from initial coords
	# then smooth up to penalty of 1000 between 5A and 15A from initial coords
	#func = 'FADE -5.0 5.0 2.0 -1000.0 1000.0'
	# FADE: lb ub d wd [ wo ] 
	constraint_text = ''
	
	RNA_residues = [ 'A', 'U', 'G', 'C']
	protein_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
		'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
		'THR', 'TRP', 'TYR', 'VAL']
	
	RNA_atoms = [ "C1'" ]
	#RNA_atoms = ["P", "C5'", "C1'", "C3'", "N1"]
	protein_atoms = [ 'CA' ]
	
	# use the RNA as ref for the protein
	# and protein as ref for the RNA
	
	RNA_ref_atom_name = ''
	prot_ref_atom_name = ''
	# first get the ref atoms
	for line in open( ref_pdb ):
		if not line.startswith('ATOM'): continue
		atom_name = line.split()[2]
		if atom_name == "C1'" and RNA_ref_atom_name == '':
			#get the RNA ref atom
			resnum = line[22:26].replace(' ','')
			chain = line[20:22].replace(' ','')
			RNA_ref_atom_name = "C1'"
			RNA_ref_resnum = resnum
			RNA_ref_chain = chain
		if atom_name == 'CA' and prot_ref_atom_name == '':
			resnum = line[22:26].replace(' ','')
			chain = line[20:22].replace(' ','')
			if chain in [ str(x) for x in range(0,10)]: continue
			prot_ref_atom_name = 'CA'
			prot_ref_resnum = resnum
			prot_ref_chain = chain
		if RNA_ref_atom_name != '' and prot_ref_atom_name != '':
			break

	for line in open( ref_pdb ):
		if not line.startswith('ATOM'): continue
		atom_name = line[12:16].replace(' ','')
		res_name = line[16:20].replace(' ','')
		is_RNA = False
		is_protein = False
		if res_name in protein_residues: is_protein = True
		if res_name in RNA_residues: is_RNA = True
		if ( not is_RNA and not is_protein) or ((atom_name not in RNA_atoms) and (atom_name not in protein_atoms)): continue
		# figure out the ref atom
		x = line[30:38].replace(' ','')
		y = line[38:46].replace(' ','')
		z = line[46:54].replace(' ','')
		resnum = line[22:26].replace(' ','')
		chain = line[20:22].replace(' ','')
		if chain in [ str(k) for k in range(0,10)]: continue
		ID = chain+resnum
		if ID in exclude_resi: continue
		if is_RNA and atom_name in RNA_atoms:
		#if atom_name == 'P':
			constraint_text += 'CoordinateConstraint %s %s%s' %(atom_name, resnum, chain)
			constraint_text += ' %s %s%s' %(prot_ref_atom_name, prot_ref_resnum, prot_ref_chain)
		elif is_protein and atom_name in protein_atoms:  # atom_name == 'CA'
			constraint_text += 'CoordinateConstraint %s %s%s' %(atom_name, resnum, chain)
			constraint_text += ' %s %s%s' %(RNA_ref_atom_name, RNA_ref_resnum, RNA_ref_chain)
		constraint_text += ' %s %s %s ' %(x,y,z)
		constraint_text += '%s\n' %(func)

	return constraint_text

def get_residues_in_range( residue_list ):
	remodel_residues = []
	all_remodel_start_res = []
	all_remodel_stop_res = []
	for resis in residue_list:
		chain = resis.split(':')[0]
		if len(resis.split('-')) < 2:
			r = resis.split(':')[1]
			remodel_residues.append( [chain, r] )
			continue
		start_res = int(resis.split(':')[1].split('-')[0])
		stop_res = int(resis.split(':')[1].split('-')[1])
		all_remodel_start_res.append( [chain, start_res] )
		all_remodel_stop_res.append( [chain, stop_res] )
		for r in range(start_res, stop_res+1):
			remodel_residues.append( [chain, r])

	return remodel_residues, all_remodel_start_res, all_remodel_stop_res

def setup_job( args ):
	# parse the pdb file
	removechain = False
	( sequences, chains, resnums) = get_sequences( args.full_struct, removechain )
	( sequences_start, chains_start, resnums_start) = get_sequences( args.start_struct, removechain )
	# figure out garbage coords (residues in the full_struct but not in the start_struct)
	garbage_coords = []
	all_full_resis = []
	all_start_resis = []
	for i,resis in enumerate(resnums_start):
		for j, r in enumerate(resis):
			chain = chains_start[i][j]
			all_start_resis.append( chain+str(r) )
	for i,resis in enumerate(resnums):
		for j, r in enumerate(resis):
			chain = chains[i][j]
			all_full_resis.append( chain+str(r) )
	garbage_coords = set(all_full_resis) - set(all_start_resis)

	# if residues are missing, fill in (skip for now)

	( coords, pdb_lines, sequence, chains2, residues2 ) = read_pdb( args.full_struct )
	# get residues within distance cutoff
	remodel_residues = []
	all_remodel_start_res = []
	all_remodel_stop_res = []
	remodel_residues, all_remodel_start_res, all_remodel_stop_res = get_residues_in_range( args.residues_to_model )
	include_in_fixed, _, _ = get_residues_in_range( args.include_in_fixed )
#	for resis in args.residues_to_model:
#		chain = resis.split(':')[0]
#		if len(resis.split('-')) < 1:
#			r = resis.split(':')[1]
#			remodel_residues.append( [chain, r] )
#			continue
#		start_res = int(resis.split(':')[1].split('-')[0])
#		stop_res = int(resis.split(':')[1].split('-')[1])
#		all_remodel_start_res.append( [chain, start_res] )
#		all_remodel_stop_res.append( [chain, stop_res] )
#		for r in range(start_res, stop_res+1):
#			remodel_residues.append( [chain, r])

	# A:1-10,A:11-12 B:1-5  --> 2 fixed chunks
	fixed_chunk_residues = [] # this is a list (each list is a diff chunk) of lists of residues
	for chunk in args.fixed_res_chunk:
		if ',' not in chunk:
			residues, _, _ = get_residues_in_range( [chunk] )
		else:
			residues, _, _ = get_residues_in_range( chunk.split(',') )
		fixed_chunk_residues.append( residues )

	CUTOFF = args.dist_cutoff
	CUTOFF2 = CUTOFF**2
	residues_near = []
	# loop through all the residues in the full structure
	# for each, loop through all the remodel residues, check if in distance cutoff
	for chain in coords.keys():
		for rsd in coords[ chain ].keys():
			if [ chain, rsd ] in remodel_residues: continue
			if chain+str(rsd) in garbage_coords: continue
			added_rsd = False
			for atom in coords[ chain ][ rsd ].keys():
				x1 = coords[ chain ][ rsd ][ atom ]
				# use C1' as a proxy
				for r_ref in remodel_residues:
					x2 = coords[ r_ref[0]][ r_ref[1]][ " C1'"]
					dist2 = 0.0
					for k in range(3):
						d = x1[k] - x2[k]
						dist2 += d*d
					if dist2 <  CUTOFF2:
						residues_near.append( [chain, rsd] )
						added_rsd = True
					if added_rsd: break # don't need to go through more atoms
				if added_rsd: break # don't need to go through more atoms

	print residues_near
	for r in include_in_fixed:
		if [r[0],int(r[1])] not in residues_near:
			print r
			residues_near.append( r )


	for line in open(args.secstruct ):
		secstruct = line.replace('\n','')
		break

	# map chain, resnum to pdb numbering from 1
	pdb_num_to_num_from_0 = {}
	index_to_pdb_num = {}
	index = 0
	full_resnum_chain_list = []
	for i,rlist in enumerate(resnums):
		for j, r in enumerate(rlist):
			pdb_num_to_num_from_0[ chains[i][j]+str(r) ] = index
			index_to_pdb_num[ index ] = [chains[i][j],r]
			index += 1
			full_resnum_chain_list.append( [chains[i][j],r] )

	remodel_indices = []
	# get the secstruct for the remodel residues
	for r in remodel_residues:
		rstr = r[0]+str(r[1])
		index = pdb_num_to_num_from_0[rstr]
		remodel_indices.append( index )

	# get the indices for the residues in the fixed chunks
	fixed_chunk_indices = []
	fixed_chunk_indices_flat = []
	for chunk in fixed_chunk_residues:
		residues_ind = []
		for r in chunk:
			rstr = r[0]+str(r[1])
			index = pdb_num_to_num_from_0[rstr]
			residues_ind.append( index )
			fixed_chunk_indices_flat.append( index )
			
		fixed_chunk_indices.append( residues_ind ) 

	# get the indices for the residues "near"
	near_indices = []
	for r in residues_near:
		rstr = r[0]+str(r[1])
		index = pdb_num_to_num_from_0[rstr]
		near_indices.append( index )

	subset_secstruct = ''
	sequence_subset = ''
	full_seq = ''
	for s in sequences:
		full_seq += s
	subset_resnum_list = []
	subset_chain_list = []
	near_only_chains = []
	near_only_resis = []
	for ind in sorted( remodel_indices + near_indices ):
		if ind in remodel_indices:
			subset_secstruct += secstruct[ind]
		else:
			near_only_chains.append( full_resnum_chain_list[ind][0] )
			near_only_resis.append( full_resnum_chain_list[ind][1] )
			subset_secstruct += '.'
		sequence_subset += full_seq[ind]
		c = full_resnum_chain_list[ind][0]
		r = full_resnum_chain_list[ind][1]
		subset_resnum_list.append( r )
		subset_chain_list.append( c )
	subset_fasta_tag = make_tag_with_dashes( subset_resnum_list, subset_chain_list )

	print subset_secstruct
	print sequence_subset

	# write secstruct
	with open('secstruct_%s.txt' %(args.job_name),'w') as f:
		f.write( subset_secstruct )
		f.write('\n')
	# write fasta
	with open('fasta_%s.txt' %(args.job_name),'w') as f:
		f.write('>%s %s\n' %(args.job_name, subset_fasta_tag))
		f.write( sequence_subset + '\n')

	# write -initial_structures file
	os.system('pdbslice.py %s -subset %s init_struct_%s_' %(args.start_struct,subset_fasta_tag,args.job_name))
	pairs_dict = dict(rna_helper.find_pairs_dict( secstruct ).items() + rna_helper.find_pairs_dict( secstruct, ['[',']'] ).items())
	print rna_helper.find_pairs_dict( secstruct, ['[',']'] ).items()
	print rna_helper.find_pairs_dict( secstruct, ['(',')'] ).items()
	print pairs_dict
#d4 = dict(d1.items() + d2.items() + d3.items())
	# get helices:
	# get continuous stretches of '('
	helices = []
	for helix_char in ['(','[','{']:
		reading_helix = False
		for i,s in enumerate(secstruct):
			if s == helix_char and not reading_helix:
				reading_helix = True
				h_start = []
			if reading_helix and s != helix_char:
				helices.append( h_start )
				reading_helix = False
			if reading_helix:
				h_start.append( i )
	split_helices = []
	for h in helices:
		stop_res = []
		for r in h:
			stop_res.append( pairs_dict[r] )
		continuous = True
		sorted_stop_res = sorted(stop_res, reverse=True)
		break_resis = []
		for i,r in enumerate(sorted_stop_res):
			if i == len( sorted_stop_res )-1: break
			if r-1 != sorted_stop_res[i+1]:
				continuous = False
				break_resis.append( i )
		if continuous:
			split_helices.append( h )
		else:
			for i,split_pos in enumerate(break_resis):
				if i == 0:
					split_helices.append( h[0:split_pos+1] )
				else:
					split_helices.append( h[break_resis[i-1]:split_pos+1] )
				if i == len( break_resis )-1:
					split_helices.append( h[split_pos+1:] )

	all_remodel_start_ind = []
	all_remodel_stop_ind = []
	for res in all_remodel_start_res:
		rstr = res[0]+str(res[1])
		index = pdb_num_to_num_from_0[rstr]
		all_remodel_start_ind.append( index )
	for res in all_remodel_stop_res:
		rstr = res[0]+str(res[1])
		index = pdb_num_to_num_from_0[rstr]
		all_remodel_stop_ind.append( index )


	# index_to_pdb_num
	additional_fixed_helix = []
	# Make helix PDBs:
	all_remodel_start_stop_ind = all_remodel_start_ind + all_remodel_stop_ind
	helix_pdb_files_to_include = []
	print "fixed_chunk_indices_flat"
	print fixed_chunk_indices_flat
	print "fixed_chunk_indices"
	print fixed_chunk_indices

	for hnum, helix in enumerate(split_helices):
		helix_resis = []
		helix_chains = []
		if helix[0] == min(remodel_indices) and pairs_dict[helix[0]] == max(remodel_indices):
			# this is the first thing getting "remodeled"
			# don't include as separate -s, instead include in the initial fixed -s structure
			for h in helix:
				additional_fixed_helix.append( h )
			continue

		if ((helix[0] in all_remodel_start_stop_ind) and (pairs_dict[helix[0]] in all_remodel_start_stop_ind)) or ((helix[-1] in all_remodel_start_stop_ind) and (pairs_dict[helix[-1]] in all_remodel_start_stop_ind)):
			for h in helix:
				additional_fixed_helix.append( h )
			continue

		if args.no_dock_connected_helix:
			if ((helix[0] in all_remodel_start_stop_ind) or (pairs_dict[helix[0]] in all_remodel_start_stop_ind)) or ((helix[-1] in all_remodel_start_stop_ind) or (pairs_dict[helix[-1]] in all_remodel_start_stop_ind)):
				for h in helix:
					additional_fixed_helix.append( h )
				continue
		# check whether it contains any fixed chunk res, then we don't want this as an input pdb
		do_not_include = False
		for hr in helix:
			if hr in fixed_chunk_indices_flat:
				print "YES"
				do_not_include = True
				break

		if do_not_include: continue
		
		for i in helix:
			# we only want the helices that are actually being remodeled
			if index_to_pdb_num[i] not in remodel_residues: continue
			helix_resis.append( index_to_pdb_num[i][1] )
			helix_resis.append( index_to_pdb_num[pairs_dict[i]][1] )
			helix_chains.append( index_to_pdb_num[i][0] )
			helix_chains.append( index_to_pdb_num[pairs_dict[i]][0] )
		if len(helix_resis) < 1: continue
		helix_tag = make_tag_with_dashes( helix_resis, helix_chains )
		#### HAVE OPTION TO NOT FIX HELIX ENDS ####
		#%s_helix_%d.pdb %(args.job_name, hnum)
		# get it from the full pdb, if it's in there
		# if not, then make it (?) - no not compatible with dock_each_chunk_per_chain, need each helix to really be positioned in the starting structure, or else it can't be a chunk
		os.system('pdbslice.py %s -subset %s TMP_HELIX_' %(args.start_struct, helix_tag))
		os.system('mv TMP_HELIX_%s %s_helix_%d.pdb' %(args.start_struct, args.job_name, hnum))
		helix_pdb_files_to_include.append( '%s_helix_%d.pdb' %(args.job_name,hnum))

	# add the fixed_chunk_res to the -s structures
	for i, chunk in enumerate(fixed_chunk_residues):
		ch_resis = []
		ch_chains = []
		for r in chunk:
			ch_resis.append( int(r[1]) )
			ch_chains.append( r[0] )
		print ch_resis
		print ch_chains
		tag = make_tag_with_dashes( ch_resis, ch_chains )
		os.system('pdbslice.py %s -subset %s TMP_fixed_chunk_' %(args.start_struct, tag))
		os.system('mv TMP_fixed_chunk_%s %s_fixed_chunk_%d.pdb' %(args.start_struct, args.job_name, i))
		# a bit of a misnomer, add to the helix_pdb_files_to_include
		helix_pdb_files_to_include.append( '%s_fixed_chunk_%d.pdb' %(args.job_name, i))

	# Add a helix to the -s structure if necessary
	n_remodel_res = len(remodel_indices)
	if len( additional_fixed_helix ) > 0:
		for ind in additional_fixed_helix:
			near_only_resis.append( full_resnum_chain_list[ind][1] )
			near_only_resis.append( full_resnum_chain_list[pairs_dict[ind]][1] )
			near_only_chains.append( full_resnum_chain_list[ind][0] )
			near_only_chains.append( full_resnum_chain_list[pairs_dict[ind]][0] )
		n_remodel_res -= 2*len(additional_fixed_helix)

	subset_near_only_tag = make_tag_with_dashes( near_only_resis, near_only_chains )
	os.system('pdbslice.py %s -subset %s near_fixed_%s_' %(args.start_struct,subset_near_only_tag,args.job_name))

	# Write coordinate constraints for all helix PDBs (and for all protein residues?)
	# make a reference PDB to use for generating constraints
	with open('ref_pdb_for_coord_csts_%s.pdb' %(args.job_name),'w') as f:
		for line in open('near_fixed_%s_%s' %(args.job_name, args.start_struct)):
			f.write(line)
		for hfil in helix_pdb_files_to_include:
			for line in open( hfil ):
				f.write(line)
	# OK and now use this to get coord csts
	if not args.no_csts:
		constraint_text = get_coord_csts( 'ref_pdb_for_coord_csts_%s.pdb' %(args.job_name))

		with open('coord_csts_%s.txt' %(args.job_name),'w') as f:
			f.write( constraint_text )

	# initial_structures file will be: init_struct_%s_%s %(args.job_name, args.start_struct)
	# write -s files
	# figure out the helix residues in the remodel section
	# write flags file
	write_flags_file( args.job_name, args.start_struct, helix_pdb_files_to_include, args.map_file, args.map_reso, n_remodel_res, args.dock_into_density, args.no_initial_structures, args.no_csts, args.extra_flags )

##################################
def write_flags_file( name, full_start_struct, 
		helix_pdb_files_to_include, 
		map_file, map_reso, 
		n_remodel_res, dock_into_density, 
		no_initial_structures, no_csts,
		extra_flags ):
	helix_string = ''
	helical_substruct_string = ''
	for h in helix_pdb_files_to_include:
		helix_string += h
		helix_string += ' '
		if "helix" in h:
			helical_substruct_string += h
			helical_substruct_string += ' '
	# get number of cycles
	cycles = int(n_remodel_res*800) # kind of made up

	with open('flags_%s' %(name), 'w') as f:
		f.write('-fasta fasta_%s.txt\n' %(name))
		f.write('-secstruct_file secstruct_%s.txt\n' %(name))
		if not no_initial_structures:
			f.write('-initial_structures init_struct_%s_%s\n' %(name,full_start_struct))
		f.write('-s near_fixed_%s_%s %s\n' %(name,full_start_struct,helix_string)) 
		f.write('-edensity:mapfile %s\n' %(map_file))
		f.write('-edensity:mapreso %s\n' %(map_reso))
		if not no_csts:
			f.write('-cst_file coord_csts_%s.txt\n' %(name))
		f.write('-new_fold_tree_initializer true\n')
		f.write('-minimize_rna true\n')
		f.write('-minimize_protein_sc true\n')
		f.write('-nstruct 500\n')
		f.write('-out:file:silent %s.out\n' %(name))
		f.write('-rna:denovo:lores_scorefxn rna/denovo/rna_lores_with_rnp_aug.wts\n') 
		f.write('-cycles %d\n' %(cycles))
		f.write('-ignore_zero_occupancy false\n')
		f.write('-output_lores_silent_file true\n')
		f.write('-rna_protein_docking true\n')
		f.write('-rnp_high_res_cycles 2\n')
		f.write('-minimize_rounds 3\n')
		f.write('-rnp_min_first true\n')
		f.write('-rnp_pack_first true\n')
		if dock_into_density:
			f.write('-dock_into_density true\n')
		else:
			f.write('-dock_into_density false\n')
		f.write('-ignore_zero_occupancy false\n')
		f.write('-convert_protein_CEN false\n')
		f.write('-FA_low_res_rnp_scoring true\n')
		f.write('-ramp_rnp_vdw true\n')
		f.write('-docking_move_size 0.3\n')
		f.write('-dock_each_chunk_per_chain true\n')
		if len( helical_substruct_string ) > 0: 
			f.write('-helical_substructs %s\n' %(helical_substruct_string))
		f.write('-mute protocols.moves.RigidBodyMover\n')
		f.write('-mute protocols.rna.denovo.movers.RNA_HelixMover\n')
		if len( extra_flags ) > 0:
			for flag in extra_flags:
				f.write( flag + ' ' )
			f.write('\n')


	# write command
	with open('command','w') as f:
		f.write('~/src/Rosetta/main_model_RNP/source/bin/rna_denovo @flags_%s' %(name))
		

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Set up sub-region RNA density modeling jobs")
	parser.add_argument('-f', '--fasta', type=str, default="", help='fasta file for the full structure')
	parser.add_argument('-d', '--dist_cutoff', type=float, default=20.0, help='distance cutoff for residues to include near remodel res')
	parser.add_argument('-ss', '--secstruct', type=str, default="", help='secstruct file for the full structure')
	parser.add_argument('-extra_flags', '--extra_flags', type=str, default="", nargs='+', help='extra flags for the rosetta run')
	parser.add_argument('-s', '--start_struct', type=str, default="", help='starting structure - must include all protein structures and at least one RNA residue')
	parser.add_argument('-full_struct', '--full_struct', type=str, default="", help='starting structure - with all the RNA with maybe garbage coords for missing RNA')
	parser.add_argument('-map', '--map_file', type=str, default="", help='map file')
	parser.add_argument('-map_reso', '--map_reso', type=str, default="", help='map resolution')
	parser.add_argument('-r', '--residues_to_model', type=str, default="", nargs='+', help='starting structure - must include all protein structures and at least one RNA residue')
	parser.add_argument('-fixed_res_chunk', '--fixed_res_chunk', type=str, default="", nargs='+', help='residues within residues_to_model that should be kept as a fixed chunk (i.e. allowed to dock, but not remodel) - specify separate chunks as e.g. A:1-10,A:11-12 B:1-5 for two chunks')
	parser.add_argument('-n', '--job_name', type=str, default="", help='name for the job')
	parser.add_argument('--dock_into_density', action='store_true', default=False, help='Dock into density - careful, your final sub-structures may no longer match up with the rest of the input structure')
	parser.add_argument('--no_initial_structures', action='store_true', default=False, help='if remodeling, do not use the starting structure as the initial coordinates')
	parser.add_argument('--no_dock_connected_helix', action='store_true', default=False, help='do not dock helices if connected by even just one residue to a fixed portion')
	parser.add_argument('--no_csts', action='store_true', default=False, help='do not use coordinate constraints to restrain positions of RNA helices/proteins')
	parser.add_argument('--include_in_fixed', type=str, default="", nargs='+',help='residues that should be in the near fixed structure')
	args = parser.parse_args()
	setup_job( args )
