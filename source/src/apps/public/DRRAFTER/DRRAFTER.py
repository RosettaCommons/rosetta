#!/usr/bin/env python
import os
import argparse
from get_sequence import get_sequences
from read_pdb import read_pdb
from make_tag import make_tag_with_dashes
import parse_tag

def get_coord_csts( ref_pdb, dist, exclude_resi=[] ):

	lb = dist * -1.0 - 5.0
	ub = dist * 1.0 + 5.0
	d = 5.0
	wd = -100.0
	wo = 100.0
	func = 'FADE %0.2f %0.2f %0.2f %0.2f %0.2f' %(lb, ub, d, wd, wo)
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

	if  RNA_ref_atom_name == '' or prot_ref_atom_name == '':
		print "ERROR: Problem generating constraints!"
		print "Your reference PDB (", ref_pdb, \
			") probably doesn't contain both protein and RNA residues."
		print "Please check your inputs!"
		print "OR:"
		print "1. Consider omitting coordinate constraints with -no_csts"
		print "2. Check out -ref_pdb_for_coord_csts"
		print "3. Check out -constrain_rigid_bodies_only (you need to have" \
			" rigid bodies that contain both protein and RNA for this!!!)"
		print "4. Create a custom constraint file and provide with -extra_flags"
		exit( 1 )

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

def find_pairs( secstruct, helix_chars=['(',')'] ):
	ss = list(secstruct) # turn into a list, so we can do replacements
	all_pairs = []  # list of dictionaries

	# find pairs
	for i in range(len(ss)):
		count = 0
		if ss[i]==helix_chars[0]:
			count += 1
			pair = {}
			pair['start'] = i
			for j in xrange(i+1,len(ss)):
				if ss[j]==helix_chars[0]: count+=1
				elif ss[j] ==helix_chars[1]: count-=1
				if count == 0: 
					pair['end'] = j
					all_pairs.append(pair)
					# get rid of those parentheses now
					ss[pair['start']] = ''
					ss[pair['end']] = ''
					break
	
	return all_pairs	

def find_pairs_dict( secstruct, helix_chars=['(',')'] ):
	lpairs = find_pairs( secstruct, helix_chars )
	pair_dict = {}
	for p in lpairs:
		pair_dict[p['start']]=p['end']
		pair_dict[p['end']]=p['start']

	return pair_dict

def parse_fasta( fasta ):
	full_chains = []
	full_resnums = []
	full_seq = []
	for line in open( fasta ):
		if line.startswith('>'):
			# get the tag from this line
			tag = ''
			for word in line.split():
				if ':' in word:
					tag += word
					tag += ' '
				
			# get the chain and residue numbers
			resnums, chains = parse_tag.parse_tag( tag )
			for c in chains:
				full_chains.append( c )
			for r in resnums:
				full_resnums.append( r )
		else:
			reading_three_letter = False
			for x in line: 
				if x == ' ': continue
				if x == '\n': continue
				if x == '[':
					reading_three_letter = True
				if reading_three_letter:
					full_seq[-1] += x
				if x == ']':
					reading_three_letter = False
				elif not reading_three_letter:
					full_seq.append( x )

	# check that the sequence is the same length as the list of chains and residue numbers
	if len(full_chains) != len(full_resnums) or \
			len(full_resnums) != len(full_seq):
		print "ERROR: Problem with the fasta file: the length of "\
			"the sequence is not the same as the length of chains/residues in tag!"
		print "Chains: ", full_chains
		print "Chain length: ", len(full_chains)
		print "Resnums: ", full_resnums
		print "Resnums length: ", len(full_resnums)
		print "Sequence: ", full_seq
		print "Sequence length: ", len(full_seq)
		exit( 1 )

	return full_chains, full_resnums, full_seq

def parse_secstruct( secstruct_file ):

	full_secstruct = []
	acceptable_ss_chars = ['.','(',')','[',']','{','}']

	for line in open( secstruct_file ):
		for x in line:
			if x not in acceptable_ss_chars:
				continue
			full_secstruct.append( x )

	return full_secstruct

def get_near_residues( struct, resnums, chains, dist_cutoff ):
	( coords, pdb_lines, sequence, chains2, residues2 ) = read_pdb( struct )

	CUTOFF = dist_cutoff
	CUTOFF2 = CUTOFF**2
	residues_near = []

	remodel_residues = []
	for x in range(len(resnums)):
		remodel_residues.append( [ chains[x], resnums[x] ] )

	for chain in coords.keys():
		for rsd in coords[ chain ].keys():
			if [ chain, rsd ] in remodel_residues: 
				residues_near.append( [chain, rsd ] )
				continue
			added_rsd = False
			for atom in coords[ chain ][ rsd ].keys():
				x1 = coords[ chain ][ rsd ][ atom ]
				# use C1' as a proxy
				for r_ref in remodel_residues:
					if " C1'" in coords[ r_ref[0]][ r_ref[1]].keys():
						x2 = coords[ r_ref[0]][ r_ref[1]][ " C1'"]
					elif " CA " in coords[ r_ref[0]][ r_ref[1]].keys(): 
						x2 = coords[ r_ref[0]][ r_ref[1]][ " CA "]
					else: continue
					dist2 = 0.0
					for k in range(3):
						d = x1[k] - x2[k]
						dist2 += d*d
					if dist2 <  CUTOFF2:
						residues_near.append( [chain, rsd] )
						added_rsd = True
					if added_rsd: break # don't need to go through more atoms
				if added_rsd: break # don't need to go through more atoms

	resnums = []
	chains = []
	for res in residues_near:
		chains.append( res[0] )
		resnums.append( res[1] )

	return resnums, chains

def setup_job( args ):

	###################################################
	# check that we have all necessary inputs
	# check here instead of letting argparse handle it
	# because these things aren't required if doing error estimation
	###################################################

	if args.map_file == "":
		print "ERROR: Please provide a -map_file!"
		exit( 1 )

	if args.map_reso == "":
		print "ERROR: Please provide -map_reso!"
		exit( 1 )

	if args.start_struct == "":
		print "ERROR: Please provide -start_struct!"
		exit( 1 )

	if args.fasta == "":
		print "ERROR: Please provide -fasta"
		exit( 1 )

	if args.secstruct == "":
		print "ERROR: Please provide -secstruct"
		exit( 1 )

	################################
	# parse the input fasta file
	################################

	full_chains, full_resnums, full_seq = parse_fasta( args.fasta )
	full_chains_resnums = []
	for x in range( len(full_chains) ):
		full_chains_resnums.append( full_chains[x] + str(full_resnums[x]) )

	full_secstruct = parse_secstruct( args.secstruct )

	# check that the secstruct is the same length as the sequence
	if len( full_secstruct ) != len( full_seq ):
		print "ERROR: Sequence and secondary structure are not the same length!"
		print "Sequence: ", full_seq
		print "Secondary structure: ", full_secstruct
		exit( 1 )


	full_secstruct_str = ''
	for x in full_secstruct:
		full_secstruct_str += x

	pairs_dict = dict(find_pairs_dict( full_secstruct_str ).items() + 
			find_pairs_dict( full_secstruct_str, ['[',']'] ).items() + 
			find_pairs_dict( full_secstruct_str, ['{','}'] ).items() )

	################################
	# figure out which residues to include in the run
	################################

	# list of indices for full_chains/full_resnums/full_seq/full_secstruct
	residues_to_include = []

	# Residues that need to be included:
	# 1. all remodel residues (should be RNA only)
	# 2. all residues that are in include_as_rigid_body_structures
	# 3. residues within distance cutoff of remodel residues 
	#    in start struct and near the -include_residues_around residues
	#    (make sure for any extra RNA residues that we include, if they 
	#    are base paired, then we also include the partner)

	##########################
	# 1. all remodel residues
	##########################
	remodel_resnums, remodel_chains = parse_tag.parse_tag( args.residues_to_model )
	remodel_chains_resnums = []
	for x in range( len(remodel_chains) ):
		remodel_chains_resnums.append( remodel_chains[x] + str( remodel_resnums[x] ) )

	for x in remodel_chains_resnums:
		residues_to_include.append( full_chains_resnums.index( x ) )

	##########################
	# 2. all residues that are in include_as_rigid_body_structures
	##########################
	rigid_body_chains_resnums = [] 
	helix_chunks = []
	for pdb in args.include_as_rigid_body_structures:
		# get the resnums and chains in the pdb
		( sequences, chains_sep, resnums_sep) = get_sequences( pdb, removechain=False )
		# each chain is its own list, put into a single list
		chains = []
		resnums = []
		for c in chains_sep:
			chains += c
		for r in resnums_sep:
			resnums += r
		
		secstruct_chunk = ''
		res_indices = []
		chains_resnums = []
		for x in range(len( chains )):
			chain_resnum = chains[x] + str( resnums[x] )
			rigid_body_chains_resnums.append( chain_resnum )
			res_index = full_chains_resnums.index( chain_resnum )
			res_indices.append( res_index )
			if res_index not in residues_to_include:
				residues_to_include.append( res_index )
			secstruct_chunk += full_secstruct[ res_index ]

		# while we're here, also figure out if this chunk is a helix
		if '.' in secstruct_chunk: continue
		is_helix_chunk = True
		for x in res_indices:
			if x not in pairs_dict.keys(): 
				is_helix_chunk = False
				break
			if pairs_dict[x] not in res_indices: 
				is_helix_chunk = False
				break
		if is_helix_chunk:
			helix_chunks.append( pdb )

	##########################
	# 3. residues within distance cutoff of remodel residues in start_struct
	##########################
	#    and near -include_residues_around
	all_include_residues_around_resnums, all_include_residues_around_chains = parse_tag.parse_tag( args.include_residues_around )
	# and then add all the remodel residues that are in the start_struct
	start_struct_seq, start_struct_chains_sep, start_struct_resnums_sep = get_sequences( args.start_struct, removechain=False )
	start_struct_chains = []
	start_struct_resnums = []
	for c in start_struct_chains_sep:
		start_struct_chains += c
	for r in start_struct_resnums_sep:
		start_struct_resnums += r

	start_struct_chains_resnums = []
	for x in range( len( start_struct_chains ) ):
		start_struct_chains_resnums.append( start_struct_chains[x] + str( start_struct_resnums[x] ) )

	for index, res in enumerate(remodel_chains_resnums):
		if res in start_struct_chains_resnums:
			all_include_residues_around_chains.append( remodel_chains[index] )
			all_include_residues_around_resnums.append( remodel_resnums[index] )

	# find the residues in the starting structure that are near these residues
	near_resnums, near_chains = get_near_residues( args.start_struct, 
			all_include_residues_around_resnums, 
			all_include_residues_around_chains, 
			args.dist_cutoff )

	near_chains_resnums = []
	for x in range( len( near_resnums ) ):
		near_chains_resnums.append( near_chains[x] + str( near_resnums[x] ) )
	
	# add the near residues to the residues_to_include
	for r in near_chains_resnums:
		res_index = full_chains_resnums.index( r )
		if res_index not in residues_to_include:
			residues_to_include.append( res_index )

	# check if all residues have a pair -- if not, include those residues as well
	extra_residues_that_pair_with_included_residues = []

	for r in residues_to_include:
		if r in pairs_dict.keys():
			# it has a partner!
			# is that partner included?
			if pairs_dict[r] not in residues_to_include:
				extra_residues_that_pair_with_included_residues.append( pairs_dict[r] )
	for r in extra_residues_that_pair_with_included_residues:
		residues_to_include.append( r )
		near_resnums.append( full_resnums[r] )
		near_chains.append( full_chains[r] )


	# also get a subset PDB from the start_struct that exludes 
	# any remodel residues -- to be included as an additional fixed chunk

	near_resnums_for_subset_pdb = []
	near_chains_for_subset_pdb = []
	for x in range( len( near_resnums ) ):
		residue = near_chains[x] + str( near_resnums[x] )
		if residue not in start_struct_chains_resnums: continue
		if residue not in remodel_chains_resnums and residue not in rigid_body_chains_resnums:
			near_resnums_for_subset_pdb.append( near_resnums[x] )
			near_chains_for_subset_pdb.append( near_chains[x] )

	created_extra_fixed_struct = False
	if len( near_resnums_for_subset_pdb ) > 0: 
		tag_for_fixed_struct = make_tag_with_dashes( near_resnums_for_subset_pdb, near_chains_for_subset_pdb )
		# TODO: modify pdbslice.py so it can be imported!
		os.system( 'pdbslice.py %s -subset %s fixed_struct_%s_' %( args.start_struct, tag_for_fixed_struct, args.job_name) )
		os.system( 'mv fixed_struct_%s_%s fixed_struct_%s.pdb' %( args.job_name, args.start_struct, args.job_name ) )
		created_extra_fixed_struct = True

	# residues_to_include
	residues_to_include_chains = []
	residues_to_include_resnums = []
	residues_to_include_sequence = ''

	for x in sorted(set(residues_to_include)):
		residues_to_include_chains.append( full_chains[x] )
		residues_to_include_resnums.append( full_resnums[x] )
		residues_to_include_sequence += full_seq[x]

	subset_to_include_tag = make_tag_with_dashes( residues_to_include_resnums, residues_to_include_chains )
	created_init_struct_pdb = False
	
	# for the initial structure pdb (any residues that aren't in the start_struct will just be omitted)
	if not args.no_initial_structures:
		os.system( 'pdbslice.py %s -subset %s init_struct_%s_' %( args.start_struct, subset_to_include_tag, args.job_name) )
		os.system( 'mv init_struct_%s_%s init_struct_%s.pdb' %(args.job_name, args.start_struct, args.job_name) )
		created_init_struct_pdb = True
	
	##########################
	# write the fasta and secstruct files for the run
	##########################
	with open('fasta_%s.txt' %(args.job_name), 'w') as f:
		f.write( '>%s %s\n' %(args.job_name, subset_to_include_tag))
		f.write( residues_to_include_sequence + '\n' )

	with open('secstruct_%s.txt' %(args.job_name), 'w') as f:
		for x in sorted(set(residues_to_include)):
			f.write( full_secstruct[x] )
		f.write('\n')

	input_structs_for_dash_s = []

	##########################
	# figure out the order of structures to be listed as -s
	# (absolute_coordinates_rigid_body_structure should be listed first)
	##########################
	if len( args.absolute_coordinates_rigid_body_structure ) > 0:
		# check that this is actually one of the rigid_body_structures
		if args.absolute_coordinates_rigid_body_structure not in args.include_as_rigid_body_structures:
			print "ERROR: -absolute_coordinates_rigid_body_structure should be one of your "\
				"include_as_rigid_body_structures!"
			print "Or if you're not including any rigid body structures,"
			print "then you can omit -absolute_coordinates_rigid_body_structure" 
			exit( 1 )
		input_structs_for_dash_s.append( args.absolute_coordinates_rigid_body_structure )
	else:
		if not created_extra_fixed_struct:
			print "Please provide -absolute_coordinates_rigid_body_structure" \
				"(one of the rigid body structures that will" \
				" set the absolute xyz for the system!)"
			exit( 1 )
		else:
			input_structs_for_dash_s.append( 'fixed_struct_%s.pdb' %( args.job_name ) )

	for pdb in args.include_as_rigid_body_structures:
		if pdb not in input_structs_for_dash_s:
			input_structs_for_dash_s.append( pdb )

	##########################
	# write coordinate constraints
	# ref pdb should contain: all include_as_rigid_body_structures residues 
	# that are present in the start_struct
	# rigid_body_chains_resnums
	##########################
	if not args.no_csts:
		constraint_text = ''
		if args.ref_pdb_for_coord_csts:
			# get a subset of this ref pdb with only the residues 
			# that we're including in the run
			os.system( 'pdbslice.py %s -subset %s ref_pdb_custom_%s_' %( args.ref_pdb_for_coord_csts, subset_to_include_tag, args.job_name) )
			constraint_text = get_coord_csts( 'ref_pdb_custom_%s_%s' %(args.job_name, args.ref_pdb_for_coord_csts ), args.cst_dist)
		elif args.constrain_rigid_bodies_only:
			ref_pdb_chains = []
			ref_pdb_resnums = []
			for x in range( len( start_struct_chains ) ):
				residue = start_struct_chains[x] + str(start_struct_resnums[x])
				if residue in rigid_body_chains_resnums:
					ref_pdb_chains.append( start_struct_chains[x] )
					ref_pdb_resnums.append( start_struct_resnums[x] )

			ref_pdb_tag = make_tag_with_dashes( ref_pdb_resnums, ref_pdb_chains )
			os.system('pdbslice.py %s -subset %s ref_pdb_%s_' %(args.start_struct,ref_pdb_tag,args.job_name))
			constraint_text = get_coord_csts( 'ref_pdb_%s_%s' %(args.job_name, args.start_struct), args.cst_dist )
		else:
			if not created_init_struct_pdb:
				print "ERROR: No init struct pdb to use for coordinate constraints!"
				print "Your options:"
				print "1. use -constrain_rigid_bodies_only: all rigid bodies that are" \
					" present in the starting structure will be constrained"
				print "2. provide a custom reference pdb for the constraints: -ref_pdb" 
				print "3. Run without coordinate constraints: provide -no_csts"
				print "4. Create a constraint file yourself (e.g. my_constraints.txt):" \
					" provide -no_csts, and also provide -extra_flags '-cst_file my_constraints.txt'"
				exit( 1 )
			constraint_text = get_coord_csts( 'init_struct_%s.pdb' %(args.job_name), args.cst_dist)
		with open('coord_csts_%s.txt' %(args.job_name), 'w') as f:
			f.write( constraint_text )

	##########################
	# write the flags file and command and print to screen
	##########################
	write_flags_file( args, helix_chunks, created_init_struct_pdb, input_structs_for_dash_s)

	print "\n"
	print "#################################################"
	print "Done setting up DRRAFTER job."
	print "* Run with 'source ./DRRAFTER_command'"
	print "* After running, the structures will be found in %s.out" %(args.job_name)
	print "* For best results, please generate at least 2000-3000 structures."
	print "* Multiple %s.out files within a directory, e.g. DRRAFTER_output/" %(args.job_name)
	print "  can be concatenated with 'easy_cat.py DRRAFTER_output/'"
	print "* PDB files can be extracted from the %s.out file with" %(args.job_name)
	print "  'extract_lowscore_decoys.py %s.out 10'" %(args.job_name)
	print "* Coordinate error in the resulting structures can be estimated"
	print "  with this script, e.g.:"
	print "  python DRRAFTER.py -final_structures %s.out.1.pdb %s.out.2.pdb [...] -estimate_error" %(args.job_name,args.job_name)
	print "#################################################"

def get_n_cycles( remodel_residues ):
	resnums, chains = parse_tag.parse_tag( remodel_residues ) 
	return min(75000, len(resnums)*800)
	

##################################
def write_flags_file( args, helix_chunks, use_init_structs, input_structs_for_dash_s):

	helix_string = ''
	for h in helix_chunks:
		helix_string += h
		helix_string += ' '

	# get number of cycles
	if args.cycles == None:
		cycles = get_n_cycles( args.residues_to_model )
	else:
		cycles = args.cycles

	dash_s_string = ''
	for s in input_structs_for_dash_s:
		dash_s_string += s
		dash_s_string += ' '

	with open('flags_%s' %(args.job_name), 'w') as f:
		f.write('-fasta fasta_%s.txt\n' %(args.job_name))
		f.write('-secstruct_file secstruct_%s.txt\n' %(args.job_name))
		if use_init_structs:
			f.write('-initial_structures init_struct_%s.pdb\n' %(args.job_name))
		f.write('-s %s\n' %(dash_s_string))
		f.write('-edensity:mapfile %s\n' %(args.map_file))
		f.write('-edensity:mapreso %s\n' %(args.map_reso))
		if not args.no_csts:
			f.write('-cst_file coord_csts_%s.txt\n' %(args.job_name))
		f.write('-new_fold_tree_initializer true\n')
		f.write('-minimize_rna true\n')
		f.write('-minimize_protein_sc true\n')
		f.write('-out:file:silent %s.out\n' %(args.job_name))
		f.write('-rna:denovo:lores_scorefxn rna/denovo/rna_lores_with_rnp_aug.wts\n') 
		f.write('-rna_protein_docking true\n')
		f.write('-rnp_min_first true\n')
		f.write('-rnp_pack_first true\n')
		if args.demo_settings:
			f.write('-cycles 500\n')
			f.write('-rnp_high_res_cycles 0\n')
			f.write('-minimize_rounds 1\n')
			f.write('-no_filters\n')
			f.write('-nstruct 10\n')
		else: 
			f.write('-cycles %d\n' %(cycles))
			f.write('-rnp_high_res_cycles 2\n')
			f.write('-minimize_rounds 3\n')
			f.write('-nstruct 500\n')
		if args.dock_into_density:
			f.write('-dock_into_density true\n')
		else:
			f.write('-dock_into_density false\n')
		f.write('-ignore_zero_occupancy false\n')
		f.write('-convert_protein_CEN false\n')
		f.write('-FA_low_res_rnp_scoring true\n')
		f.write('-ramp_rnp_vdw true\n')
		f.write('-docking_move_size %s\n' %(args.docking_move_size))
		if args.no_dock_each_rigid_body_separately:
			f.write('-dock_each_chunk_per_chain false\n')
		else:
			f.write('-dock_each_chunk_per_chain true\n')
		if len( helix_string ) > 0: 
			f.write('-helical_substructs %s\n' %(helix_string))
		f.write('-mute protocols.moves.RigidBodyMover\n')
		f.write('-mute protocols.rna.denovo.movers.RNA_HelixMover\n')
		f.write('-use_legacy_job_distributor true\n')
		if len( args.extra_flags ) > 0:
			for flag in args.extra_flags:
				f.write( '-' + flag + ' ' )
			f.write('\n')

	# write command
	rosetta_dir = args.rosetta_directory
	if len(rosetta_dir) > 0 and not rosetta_dir.endswith('/'):
		rosetta_dir += '/'
	with open('DRRAFTER_command','w') as f:
		f.write('%srna_denovo @flags_%s\n' %(rosetta_dir, args.job_name))

def estimate_error( args ):
	# check that final structures were provided
	if len( args.final_structures ) < 1: 
		print "ERROR: Please provide the final DRRAFTER models as -final_structures"
		exit( 1 )
	elif len( args.final_structures ) < 10:
		print "WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
		print "For best results, provide at least 10 structures for error estimation."
		print "WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
	
	struct_str = ''
	for s in args.final_structures:
		struct_str += s
		struct_str += ' '
	rosetta_dir = args.rosetta_directory
	if len(rosetta_dir) > 0 and not rosetta_dir.endswith('/'):
		rosetta_dir += '/'
	command = '%sdrrafter_error_estimation -s %s -mute core' %(rosetta_dir, struct_str )
	print "\nRUNNING ", command, '\n'
	os.system( '%sdrrafter_error_estimation -s %s -mute core' %(rosetta_dir, struct_str ))


if __name__ == '__main__':

	parser=argparse.ArgumentParser( description="DRRAFTER: De novo RNP modeling in Real-space through Assembly of Fragments Together with Electron density in Rosetta" )
	parser.add_argument('-fasta', type=str, default="",
		help='fasta file for the full structure')
	parser.add_argument('-secstruct', type=str, default="",
		help='secstruct file for the full structure')
	parser.add_argument('-rosetta_directory', '--rosetta_directory',
		type=str, default="", help='Path to Rosetta executables')
	parser.add_argument('-start_struct', type=str, default="", 
		help='starting structure - must include all protein structures and at least one RNA residue')
	parser.add_argument('-residues_to_model', '-r', type=str, 
		default="", nargs='+', help='starting structure - must '
		'include all protein structures and at least one RNA residue')
	parser.add_argument('-map_file', type=str, default="", help='map file')
	parser.add_argument('-map_reso', type=str, default="", help='map resolution')
	parser.add_argument('-include_residues_around', type=str, default="", 
		nargs='+', help='residues in start_struct near the region to '
		'be modeled. Residues within dist_cutoff of these residues '
		'will be included in the run.')
	parser.add_argument('-dist_cutoff', type=float, default=20.0,
		 help='distance cutoff for residues to include near remodel res')
	parser.add_argument('-extra_flags', type=str, 
		default="", nargs='+', help='extra flags for the rosetta run')
	parser.add_argument('-docking_move_size', 
		type=str, default="0.5", help='between 0.0 and 1.0, '
		'1.0 corresponds to the most aggressive docking, 0.0 to the least aggressive')
	parser.add_argument('-include_as_rigid_body_structures', type=str, 
		default="", nargs='+', help='structures that should be '
		'treated as rigid bodies during modeling (e.g. provide '
		'separate files for each of the helices that were fit into the density)')
	parser.add_argument('-absolute_coordinates_rigid_body_structure', type=str, 
		default="", help='One of the structures that will be '
		'included as a rigid body that will set the absolute '
		'xyz coordinates for the run (this structure must be '
		'fit into the density map)')
	parser.add_argument('-job_name', '-n', type=str, default="", help='name for the job')
	parser.add_argument('-dock_into_density', action='store_true', 
		default=False, help='Dock into density - careful, '
		'your final sub-structures may no longer match up '
		'with the rest of the input structure')
	parser.add_argument('-no_dock_each_rigid_body_separately', 
		action='store_true', default=False, 
		help='Determines how the kinematics are set up')
	parser.add_argument('-no_initial_structures', action='store_true', 
		default=False, help='if remodeling, do not use '
		'the starting structure as the initial coordinates')
	parser.add_argument('-constrain_rigid_bodies_only', action='store_true', 
		default=False, help='for coordinate constraints, only '
		'constrain residues in the included rigid bodies')
	parser.add_argument('-no_csts', action='store_true', default=False, 
		help='do not use coordinate constraints to restrain positions of RNA helices/proteins')
	parser.add_argument('-cst_dist', default = 10.0, type=float, 
		help='Distance tolerance for constraints')
	parser.add_argument('-ref_pdb_for_coord_csts',
		type=str, default="", help='reference pdb for coordinate '
		'constraints (must include at least one RNA and one protein residue)')
	parser.add_argument('-cycles', type=int, help='Number of monte carlo cycles')
	parser.add_argument('-demo_settings', action='store_true', default=False, 
		help='Settings to make the demo run quickly - do NOT use in actual runs')
	parser.add_argument('-estimate_error', action='store_true', 
		default=False, help='Estimate error in final models')
	parser.add_argument('-final_structures', type=str, default="", 
		nargs='+', help='Final DRRAFTER models (e.g. top 10 scoring), used for error estimation')
	args = parser.parse_args()

	if args.estimate_error:
		estimate_error( args )
	else:
		setup_job( args )
