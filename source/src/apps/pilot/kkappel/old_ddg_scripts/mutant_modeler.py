from sys import exit
from read_pdb import read_pdb
from collections import OrderedDict
import glob
import os
from itertools import islice
from itertools import permutations
from parse_tag import parse_tag
import get_surrounding_res
from copy import deepcopy
from make_tag import make_tag_with_dashes

print "\n#################################################"
print "WARNING: Please use updated versions of these scripts in apps/public/rnp_ddg/"
print "#################################################\n"

protein_residues = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 
			'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR' ]
protein_residues_1letter = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 
			'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
RNA_residues = ['  A', '  U', '  C', '  G']
RNA_residues_1letter = [ 'a', 'u', 'c', 'g' ]
pdb_to_one_letter = { '  A': 'a',
			'  U': 'u',
			'  C': 'c',
			'  G': 'g',
			'ALA': 'A',
			'CYS': 'C',
			'ASP': 'D',
			'GLU': 'E',
			'PHE': 'F',
			'GLY': 'G',
			'HIS': 'H',
			'ILE': 'I',
			'LYS': 'K',
			'LEU': 'L',
			'MET': 'M',
			'ASN': 'N',
			'PRO': 'P',
			'GLN': 'Q',
			'ARG': 'R',
			'SER': 'S',
			'THR': 'T',
			'VAL': 'V',
			'TRP': 'W',
			'TYR': 'Y' }

def figure_out_WT_seq( wt_struct ):
	(coords, pdb_lines, sequence, chains, residues) = read_pdb( wt_struct )
	# preserve the chain order from the wt_struct
	chain_letters = list(OrderedDict.fromkeys(chains))
	# figure out which of these is protein and which is RNA
	protein_chains = []
	RNA_chains = []
	for c in chain_letters:
		for r in sequence[c].values():
			if r in protein_residues:
				protein_chains.append(c)
			elif r in RNA_residues:
				RNA_chains.append(c)
	# Make sure that there aren't chains that contain both RNA and protein residues
	if len( set(RNA_chains) & set(protein_chains) ) > 0:
		raise ValueError('Chain %s contains both RNA and protein residues' %( set(RNA_chains) & set(protein_chains)))

	# Write out a new pdb, all protein (sorted chains), then all RNA (sorted chains)
	# Make sure the residues are sorted as well
	full_residue_list = []
	full_chain_list = []
	with open(wt_struct.replace('.pdb','')+'_reorder.pdb', 'w') as newpdb:
		sorted_list_of_chains = sorted(set(protein_chains))
		for c in sorted(set(RNA_chains)):
			sorted_list_of_chains.append( c )
		for c in sorted_list_of_chains:
			residues = sorted(pdb_lines[c].keys())
			for resi in residues:
				full_residue_list.append(resi)
				full_chain_list.append(c)
				for atom in pdb_lines[c][resi].keys():
					newpdb.write( pdb_lines[c][resi][atom] + '\n')

	# figure out the wildtype sequence
	wt_seq_full = ''
	wt_RNA_seq = ''
	for c in sorted_list_of_chains:
	#for c in chain_letters:
		# Assuming that we want the residues to be in order of residue number
		for r in sorted(sequence[c].keys()):
			wt_seq_full+= pdb_to_one_letter[ sequence[c][r] ]
			if c in set(RNA_chains):
				wt_RNA_seq+= pdb_to_one_letter[ sequence[c][r] ]

	# rna nums (useful for later reference/renumbering)
	rna_nums = ''
	for c in sorted(set(RNA_chains)):
		rna_nums += c + ':' + str(min(sequence[c].keys())) + '-'
		rna_nums += str(max(sequence[c].keys())) + ' '
			
	wt_pdb_tag = make_tag_with_dashes( full_residue_list, full_chain_list )

	return set(RNA_chains), set(protein_chains), wt_seq_full, wt_RNA_seq, rna_nums, wt_pdb_tag


# Takes a sequence and returns the ss
# I'm sure there's a better way of doing this...
def get_ss( seq, T=37 ):
	# make a temporary sequence file
	ftmp = open('tmp_seq.txt', 'w')
	ftmp.write( seq +'\n')
	ftmp.close()
	os.system('RNAfold -T %d <tmp_seq.txt> tmp_ss.out' %(T))
	ss = ''
	with open( 'tmp_ss.out', 'r') as ssfile:
		for line in islice( ssfile, 1, 2):
			ss = line.split()[0]
			
	# remove the temporary sequence and output files
	os.system('rm tmp_seq.txt')
	os.system('rm tmp_ss.out')

	return ss


# returns a list of dictionaries, with start and end values for all base pairs
def find_pairs( secstruct ):
	ss = list(secstruct) # turn into a list, so we can do replacements
	all_pairs = []  # list of dictionaries

	# find pairs
	for i in range(len(ss)):
		count = 0
		if ss[i]=='(':
			count += 1
			pair = {}
			pair['start'] = i
			for j in xrange(i+1,len(ss)):
				if ss[j]=='(': count+=1
				elif ss[j] ==')': count-=1
				if count == 0: 
					pair['end'] = j
					all_pairs.append(pair)
					# get rid of those parentheses now
					ss[pair['start']] = ''
					ss[pair['end']] = ''
					break
	
	return all_pairs

# Returns a list of indices that index input bps
def list_NC_BPs( bps, #list of dicts 
		seq ):

	list_nc_bps = []
	# loop through the list of base pairs
	for i in range(len(bps)):
		if ((seq[bps[i]['start']] == 'a') and (seq[bps[i]['end']] == 'u')
			or (seq[bps[i]['start']] == 'u') and (seq[bps[i]['end']] == 'a')
			or (seq[bps[i]['start']] == 'g') and (seq[bps[i]['end']] == 'c')
			or (seq[bps[i]['start']] == 'c') and (seq[bps[i]['end']] == 'g')):
			continue
		list_nc_bps.append(i)

	# Return a list of indices of nc base pairs, can be used to index bps
	return list_nc_bps 

# Like list_NC_BPs, returns a list of indices that index input bps
def list_mutated_BPs( bps, seq ):
	list_mut_bps = []
	# loop through the list of base pairs
	for i in range(len(bps)):
		# check if base pair sequence is different from WT
		if ((seq[bps[i]['start']] != wt_seq[bps[i]['start']]) or 
			(seq[bps[i]['end']] != wt_seq[bps[i]['end']])):
			list_mut_bps.append(i)

	return list_mut_bps

class MutantModeler:
	"""
	Class in charge of making mutants of a given WT structure
	Generates commands that will make mutant structures 
	starting from a 'WT' structure
	Uses either the low-res, med-res, or high-res protocols
	"""
	def __init__( self, method, wt_struct ):
		self.wt_struct = wt_struct.replace('.pdb','') + '_reorder.pdb' #keep the name of the wt struct file
		self.method = method
		self.run_local = True
		self.mutant_seqs = []
		self.high_res_WT_rebuild = [] # a list of lists of mutated loop residues
		self.med_res_WT_rebuild = [] # a list of lists of non-canonical rebuild positions
		self.compare_to_full_WT = False # Will be reset during first call to add_mutant
		# Figure out the WT sequence and which chains are RNA and protein
		# This will also create a new reordered version of the pdb 
		# all protein will come first, then all RNA
		self.RNA_chains, self.protein_chains, self.wt_seq_full, self.wt_RNA_seq, self.RNA_nums, self.wt_pdb_tag = figure_out_WT_seq( wt_struct )
		# Figure out the WT RNA secondary structure, this can be overwritten later
		self.wt_RNA_secstruct = get_ss( self.wt_RNA_seq )
		self.mutant_index = 0
		self.cmd_opts = ''
		self.sfxn = 'P_overlap_reR-hbond_sc_MOD'
		self.protein_pack_reps = 10
		self.rosetta_prefix = ''
		self.wt_base_pairs = find_pairs( self.wt_RNA_secstruct )
		self.map_RNA_index_to_full = {}
		RNA_index = 0
		for i in range(len(self.wt_seq_full)):
			if self.wt_seq_full[i] in RNA_residues_1letter:
				self.map_RNA_index_to_full[RNA_index]=i
				RNA_index+=1
		# the master subset dictionary, currently only used for high-res
		self.mutant_subsets = []

	def set_wt_RNA_secstruct( self, secstruct ):
		self.wt_RNA_secstruct = secstruct
		# Update the base pairs using this secondary structure
		self.wt_base_pairs = find_pairs( self.wt_RNA_secstruct )

	def set_run_locally( self, value ):
		self.run_local = value

	def get_mut_positions( self, seq ):
		# compare to the WT seq
		mut_pos = []
		if self.compare_to_full_WT:
			compare_seq = self.wt_seq_full
		else: 
			compare_seq = self.wt_RNA_seq

		for i in range(0,len(seq)):
			if ( seq[i] != compare_seq[i] ):
				mut_pos.append( i )
		return mut_pos

	def add_mutant( self, seq ):
		# First mutant added, figure out if we're dealing with full sequences or just RNA sequences
		# For now this method can only handle mutations to the same length (no additions or deletions)
		if len( self.mutant_seqs ) == 0:
			if len(seq) == len( self.wt_seq_full ):
				self.compare_to_full_WT = True
				# The wildtype sequence should be have index 0
				self.mutant_seqs.append( [self.wt_seq_full, []] )
			else:
				if (len(seq) != len( self.wt_RNA_seq )):
					raise ValueError('The mutant sequence needs to be the same length as the WT  seq!')
				self.compare_to_full_WT = False
				# The wildtype sequence should be have index 0
				self.mutant_seqs.append( [self.wt_RNA_seq, []] )

		# figure out mutated positions (returns a list)
		mutated_positions = self.get_mut_positions( seq )
		self.mutant_seqs.append( [seq, mutated_positions] )
		self.mutant_index+=1
		return self.mutant_index
	
	# Define extra command line options
	def set_cmd_opts( self, cmd_opts ):
		self.cmd_opts = cmd_opts

	def set_sfxn( self, sfxn ):
		self.sfxn = sfxn

	def set_protein_pack_reps( self, protein_pack_reps ):
		self.protein_pack_reps = protein_pack_reps

	def set_rosetta_prefix( self, prefix ):
		self.rosetta_prefix = prefix

	def write_mutfile( self, mut_indices, mutant_seq, input_file, mfil_name ):
		# Write the mutfile
		with open( mfil_name, 'w') as mfil:
			total_muts = len(mut_indices)
			mfil.write('total %d \n' %(total_muts))
			mfil.write('%d\n' %(total_muts))
			for i in mut_indices:
				if self.compare_to_full_WT:
					mfil.write('%s %d %s\n' %(self.wt_seq_full[i], i+1, mutant_seq[i]))
				else:
					mfil.write('%s %d %s\n' %(self.wt_RNA_seq[i], self.map_RNA_index_to_full[i]+1, mutant_seq[i]))
			mfil.write('input_file %s\n' %(input_file))
		return

	def figure_out_protein_pack_reps( self, index ):
		# If this is the WT sequence, then do protein packing
		if index == 0:
			return self.protein_pack_reps
		# Figure out if there are protein mutations, if yes, then do protein packing
		elif self.compare_to_full_WT:
			for mut_index in self.mutant_seqs[index][1]:
				if mut_index < (len(self.wt_seq_full) - len(self.wt_RNA_seq)):
					# Then this is a protein mutation
					return self.protein_pack_reps
		# if no protein mutations and not the WT seq, then return 0 protein pack reps
		return 0
			
	
	def get_low_res_command_lines( self, index ):
		# Write the mutfile
		self.write_mutfile( self.mutant_seqs[index][1], self.mutant_seqs[index][0], 'none', 'mutfile' )

		# Now write the mutate_and_score_RNP command line
		cmd = ''
		cmd += '%smutate_and_score_RNP -ignore_zero_occupancy false ' %(self.rosetta_prefix)
		cmd += '-s %s -mutfile mutfile ' %( self.wt_struct )
		cmd += '-out_prefix ddG_ '
		pack_reps = self.figure_out_protein_pack_reps( index )
		cmd += '-protein_pack_reps %d ' %(pack_reps)
		cmd += '%s > get_complex_scores.log\n' %( self.cmd_opts )
		return cmd

	# Non-canonical base pairs will be modeled with farfar
	def get_med_res_command_lines( self, index ):
		model_cmd_files = []
		# Find the non-canonical base pairs
		if self.compare_to_full_WT:
			RNA_seq = [self.mutant_seqs[index][0][i] for i in sorted(self.map_RNA_index_to_full.values())]
			RNA_seq = ''.join(RNA_seq)
		else:
			RNA_seq = self.mutant_seqs[index][0]
		nc_base_pairs = list_NC_BPs( self.wt_base_pairs, RNA_seq )

		# If there aren't any non-canonical base pairs, then just call
		# the low-res method
		if ( len(nc_base_pairs) == 0 ): # no non-canonicals
			final_ddg_commands = self.get_low_res_command_lines( index )
			return '', model_cmd_files, final_ddg_commands, []
		# Do this part only if there are non-canonical base pairs
		# First need to thread any canonical mutations
		list_nc_pos = [] # List of indices of NC BP positions
		for i in nc_base_pairs:
			list_nc_pos.append( self.wt_base_pairs[i]['start'] )
			list_nc_pos.append( self.wt_base_pairs[i]['end'] )
		surround_pos = set() # List of indices of BPs surrounding NC BPs
		# TODO: probably want to do this differently if there are big loops btwn BPs
		for bp in nc_base_pairs:
			if bp > 0: # Add the base pair below (as long as it's not the 1st BP)
				surround_pos.add( self.wt_base_pairs[bp-1]['start'] )
				surround_pos.add( self.wt_base_pairs[bp-1]['end'] )
			if bp < len(self.wt_base_pairs)-1: 
				# Add the base pair above (as long as it's not the last base pair)
				surround_pos.add( self.wt_base_pairs[bp+1]['start'] )
				surround_pos.add( self.wt_base_pairs[bp+1]['end'] )
		# get rid of "surrounding positions" that are actually NC BPs
		surround_pos-=set(list_nc_pos) 
		# Find mutated positions that aren't NC BPs or surrounding
		other_mut_pos = [] #should be in full indexing if compare_to_full_WT, otherwise RNA indexing
		for i in self.mutant_seqs[index][1]:
			if self.compare_to_full_WT:
				# Then we need to map the full index back to the RNA index
				# b/c list_nc_pos and surround_pos and in RNA indexing, need to
				# compare to these lists below
				m_ind = i - (len(self.wt_seq_full) - len(self.wt_RNA_seq))
			else: m_ind = i
			if (m_ind in list_nc_pos) or (m_ind in surround_pos): continue
			other_mut_pos.append(i)
		# Write the mutfile to make these other mutations
		setup_cmd = ''
		if len(other_mut_pos) > 0 : # only if there are other mutations
			self.write_mutfile( other_mut_pos, self.mutant_seqs[index][0], 'none', 'mutfile_initial' )

			# Write the setup commands
			# First mutate_and_score_RNP run to make the initial mutations
			setup_cmd += '%smutate_and_score_RNP -ignore_zero_occupancy false ' %(self.rosetta_prefix)
			setup_cmd += '-s %s -mutfile mutfile_initial ' %( self.wt_struct )
			setup_cmd += '-out_prefix initial_muts_ '
			setup_cmd += '%s > make_initial_muts.log\n' %( self.cmd_opts )

			# Next the RNA chain needs to be extracted from the resulting PDB
			setup_cmd += 'cp complex_lowest_scoring_initial_muts_*.pdb complex_start_struct.pdb\n'
			#setup_cmd += 'cp initial_muts*_bound.pdb complex_start_struct.pdb\n'
			setup_cmd += 'extract_chain.py complex_lowest_scoring_initial_muts_*.pdb %s\n' %(' '.join(self.RNA_chains))

			# Rename the resulting PDB for clarity
			setup_cmd += 'mv complex_lowest_scoring_initial_muts_*%s.pdb RNA_bound_conf.pdb\n' %(''.join(self.RNA_chains))

		else: # if there aren't other mutations
			# Don't need to run mutate_and_score_RNP first
			# Need to get the RNA bound structure from the starting structure
			setup_cmd += 'cp %s complex_start_struct.pdb\n' %(self.wt_struct)
			setup_cmd += 'extract_chain.py complex_start_struct.pdb %s\n' %(' '.join(self.RNA_chains))
			setup_cmd += 'mv complex_start_struct%s.pdb RNA_bound_conf.pdb\n' %(''.join(self.RNA_chains ) )

		# Renumber the RNA_bound_conf for simplicity (will renumber again later)
		setup_cmd += 'renumber_pdb_in_place.py RNA_bound_conf.pdb A:1-%d\n' %(len(RNA_seq))
		setup_cmd += 'renumber_pdb_in_place.py RNA_bound_conf.pdb A:1-%d\n' %(len(RNA_seq))
		# Now excise the residues that should get rebuilt with farfar
		excise_str = ''
		excise_list = []
		excise_pos = surround_pos
		for i in list_nc_pos: 
			excise_pos.add( i )
		for i in excise_pos:
			excise_str+=('A'+str(i+1)+' ')
			excise_list.append(i)

		# If any of the first 3 BPs is being remodeled, rebuild it and all the BPs below it
		# TODO: is this logic sufficiently general??
		# For MS2, this definitely makes sense b/c these bps are not really interacting
		# with the protein, but it might make less sense if these residues are interacting w/ protein
		end_bp_res = [2,1] # don't need to include first basepair b/c nothing is "below" it
		cut_below_point = 0
		for i in end_bp_res: # checking for nc muts in last 3 basepairs
			if self.wt_base_pairs[i]['start'] in list_nc_pos:
				#this means that we're remodeling one of the last 3 non-canonical base pairs
				cut_below_point = i
				break # checking in highest to lowest order, so this is the top of where to cut
		if cut_below_point > 0: 
			for i in range( cut_below_point ):
				excise_str+=('A'+str(self.wt_base_pairs[i]['start']+1)+' ')
				excise_list.append(self.wt_base_pairs[i]['start'])
				excise_str+=('A'+str(self.wt_base_pairs[i]['end']+1)+' ')
				excise_list.append(self.wt_base_pairs[i]['end'])
		########## End figuring out end BPs to remodel ###############

		# Add the excision command to the setup commands
		setup_cmd += 'pdbslice.py RNA_bound_conf.pdb -excise %s scaff_\n' %(excise_str)
		setup_cmd += 'pdbslice.py RNA_bound_conf.pdb -excise %s scaff_\n' %(excise_str)

		# Write a secondary structure string for farfar, can't have non-canonical BPs
		seq_ss_list = list( self.wt_RNA_secstruct )
		for i in list_nc_pos:
			seq_ss_list[i] = '.'
		new_ss_str = ''.join(seq_ss_list)
		with open('rna_seq.fasta', 'w') as fasta:
			fasta.write( '> RNA_bound_conf.pdb A:1-%d\n' %(len(RNA_seq)))
			fasta.write( '%s\n' %(RNA_seq))

		# Make the rna_denovo input files
		# TODO: HACK FOR TESTING!!!
		# TODO: Set the number of cycles and the nstruct based on number of excised residues!!
		#setup_cmd += 'rna_denovo_setup.py -nstruct 20 -cycles 10 -fasta rna_seq.fasta -secstruct "%s" -fixed_stems ' %(new_ss_str)
		setup_cmd += 'rna_denovo_setup.py -nstruct 200 -fasta rna_seq.fasta -secstruct "%s" -fixed_stems ' %(new_ss_str)
		setup_cmd += '-chemical:enlarge_H_lj -s scaff_RNA_bound_conf.pdb -tag RNA_seq '
		setup_cmd += '-score:weights %s\n' %(self.sfxn)

		# Once this is run, it creates 'README_FARFAR'
		# TODO: get rid of this hack to change rosetta directory
		setup_cmd += 'sed -i "s/main/main_RNP/g" README_FARFAR\n'
		if self.run_local:
			setup_cmd += 'source README_FARFAR\n'
		#else:
		# TODO: option to change how long this runs?
			#setup_cmd += 'rosetta_submit.py README_FARFAR farfar_out 1 24\n'
			#setup_cmd += 'source ./qsubMINI\n'
		model_cmd_files.append('%d/README_FARFAR' %(index))

		# Write a file that lists the excised residues to be read later by mutate_and_score_RNP
		excise_list_index_full = []
		# Get these indices correct so that we don't need to use -prot_offset_num and -RNA_offset_num
		for i in excise_list:
			# Add 1 b/c residue numbering is 1 based in rosetta
			excise_list_index_full.append(str(self.map_RNA_index_to_full[i]+1))
		with open('REBUILD_RESIDUES.txt', 'w') as resfile:
			resfile.write(' '.join(excise_list_index_full)+'\n')

		# Done with the setup commands (this will run the farfar jobs
		# Now create the finalize commands, these need to cleanup results from setup
		# and submit the final ddG jobs
		
		# Write a mutfile for the final ddG job (there are no more mutations, there will just be
		# minimization, repacking, and scoring)
		with open( 'mutfile', 'w') as mfil:
			mfil.write('total 0\n0\n')
			mfil.write('input_file top20_complexes.out \n')

		# Get rid of any previous structures in this directory
		final_cmd = 'rm top20_complexes.out\n'
		final_cmd += 'rm complex.out.*.pdb\n'

		# Get the results from the farfar run
		#final_cmd += 'easy_cat.py farfar_out\n'
		final_cmd += 'extract_lowscore_decoys.py RNA_seq.out 20 \n'

		# Concatenate the new RNA structure with the protein structure
		final_cmd += 'extract_chain.py complex_start_struct.pdb %s\n' %(' '.join(self.protein_chains))
		for i in range(1,21):
			# First renumber the RNA back
			final_cmd += 'renumber_pdb_in_place.py RNA_seq.out.%d.pdb %s\n' %(i, self.RNA_nums)
			final_cmd += 'renumber_pdb_in_place.py RNA_seq.out.%d.pdb %s\n' %(i, self.RNA_nums)
			# Then concatenate
			final_cmd += 'cat complex_start_struct%s.pdb RNA_seq.out.%d.pdb > complex.out.%d.pdb\n' %(''.join(self.protein_chains),i,i)
		# Combine all the pdbs into a single silent file
		final_cmd += 'combine_silent -s complex.out.*.pdb -out:file:silent_struct_type binary '
		final_cmd += '-out:file:silent top20_complexes.out\n'

		# Now do the final mutate_and_score_RNP run to get the complex scores
		final_cmd += '%smutate_and_score_RNP -ignore_zero_occupancy false ' %(self.rosetta_prefix)
		final_cmd += '-s %s -mutfile mutfile ' %( self.wt_struct )
		final_cmd += '-out_prefix ddG_ '
		pack_reps = self.figure_out_protein_pack_reps( index )
		final_cmd += '-protein_pack_reps %d ' %(pack_reps)
		final_cmd += '%s > get_complex_scores.log\n' %( self.cmd_opts )

		################ FIGURE OUT WT REBUILD STUFF ######################

		wt_index = -1
		# Figure out if the appropriate WT structure has already been built
		for wt in range(len(self.med_res_WT_rebuild)):
			if set(list_nc_pos) == set(self.med_res_WT_rebuild[wt]):
				wt_index = wt
				break
		# If it hasn't been built, then build it here
		if wt_index == -1: # hasn't already been built
			mut_dir = os.getcwd()
			os.chdir('..')
			base_dir = os.getcwd()
			if not os.path.exists('wt_%d' %(len(self.med_res_WT_rebuild)) ):
				os.mkdir('wt_%d' %(len(self.med_res_WT_rebuild)) )
			os.chdir( 'wt_%d' %(len(self.med_res_WT_rebuild)) )
			# Copy the wt start struct here as the "native"
			os.system('cp %s complex_start_struct.pdb' %(self.wt_struct))
			# Extract the RNA chains
			os.system('extract_chain.py complex_start_struct.pdb %s' %(' '.join(self.RNA_chains)))
			# Rename for simplicity 
			os.system('mv complex_start_struct%s.pdb RNA_bound_conf.pdb' %(''.join(self.RNA_chains ) ))
			# Renumber for simplicity
			os.system('renumber_pdb_in_place.py RNA_bound_conf.pdb A:1-%d' %(len(self.wt_RNA_seq)))
			os.system('pdbslice.py RNA_bound_conf.pdb -excise %s scaff_' %(excise_str))

			with open('rna_seq.fasta', 'w') as fasta:
				fasta.write( '> RNA_bound_conf.pdb A:1-%d\n' %(len(self.wt_RNA_seq)))
				fasta.write( '%s\n' %(self.wt_RNA_seq))

			wt_rna_denovo_setup_cmd = ''
			wt_rna_denovo_setup_cmd += 'rna_denovo_setup.py -nstruct 200 -fasta rna_seq.fasta -secstruct "%s" -fixed_stems ' %(new_ss_str)
			wt_rna_denovo_setup_cmd += '-chemical:enlarge_H_lj -s scaff_RNA_bound_conf.pdb -tag RNA_seq '
			wt_rna_denovo_setup_cmd += '-score:weights %s' %(self.sfxn)
			os.system( wt_rna_denovo_setup_cmd )
			# This makes a README_FARFAR file
			os.system('sed -i "s/main/main_RNP/g" README_FARFAR')
			if self.run_local:
				setup_cmd += 'cd wt_%d\n' %(len(self.med_res_WT_rebuild))
				setup_cmd += 'source README_FARFAR\n'
			model_cmd_files.append('wt_%d/README_FARFAR' %(len(self.med_res_WT_rebuild)))
			with open('REBUILD_RESIDUES.txt', 'w') as resfile:
				resfile.write(' '.join(excise_list_index_full)+'\n')
			# Done with the "setup" (nothing gets added if not running locally b/c it just gets run in this code)

			# Write a mutfile for the final ddG job (there are no more mutations, there will just be
			# minimization, repacking, and scoring)
			with open( 'mutfile', 'w') as mfil:
				mfil.write('total 0\n0\n')
				mfil.write('input_file top20_complexes.out \n')

			# Get rid of any previous structures in this directory
			final_cmd += 'cd wt_%d\n' %(len(self.med_res_WT_rebuild))
			final_cmd += 'rm top20_complexes.out\n'
			final_cmd += 'rm complex.out.*.pdb\n'

			# Get the results from the farfar run
			#final_cmd += 'easy_cat.py farfar_out\n'
			final_cmd += 'extract_lowscore_decoys.py RNA_seq.out 20 \n'

			# Concatenate the new RNA structure with the protein structure
			final_cmd += 'extract_chain.py complex_start_struct.pdb %s\n' %(' '.join(self.protein_chains))
			for i in range(1,21):
				# First renumber the RNA back
				final_cmd += 'renumber_pdb_in_place.py RNA_seq.out.%d.pdb %s\n' %(i, self.RNA_nums)
				final_cmd += 'renumber_pdb_in_place.py RNA_seq.out.%d.pdb %s\n' %(i, self.RNA_nums)
				# Then concatenate
				final_cmd += 'cat complex_start_struct%s.pdb RNA_seq.out.%d.pdb > complex.out.%d.pdb\n' %(''.join(self.protein_chains),i,i)
			# Combine all the pdbs into a single silent file
			final_cmd += 'combine_silent -s complex.out.*.pdb -out:file:silent_struct_type binary '
			final_cmd += '-out:file:silent top20_complexes.out\n'

			# Now do the final mutate_and_score_RNP run to get the complex scores
			final_cmd += '%smutate_and_score_RNP -ignore_zero_occupancy false ' %(self.rosetta_prefix)
			final_cmd += '-s %s -mutfile mutfile ' %( self.wt_struct )
			final_cmd += '-out_prefix ddG_ '
			final_cmd += '-protein_pack_reps %d ' %(pack_reps)
			final_cmd += '%s > get_complex_scores.log\n' %( self.cmd_opts )
			
			os.chdir( mut_dir )
			wt_index = len(self.med_res_WT_rebuild)
			self.med_res_WT_rebuild.append(list_nc_pos)

		# write a wt_comparison.txt file in the dir that says which wt struct should be
		# used for comparison
		with open( 'wt_comparison.txt', 'w') as cfil:
			cfil.write( 'wt_' + str(wt_index) + '\n')

		############END WT REBUILD ##############################################

		return setup_cmd, model_cmd_files, final_cmd, []
			

	def get_high_res_command_lines( self, index ):
		model_cmd_files_with_dir = []
		# First figure out if there are any mutations in base-paired regions
		# these mutations should just be threaded
		# Mutation positions self.mutant_seqs[ index ][1]
		# wt secstruct: self.wt_RNA_secstruct
		if self.compare_to_full_WT:
			RNA_seq = [self.mutant_seqs[index][0][i] for i in sorted(self.map_RNA_index_to_full.values())]
			RNA_seq = ''.join(RNA_seq)
		else:
			RNA_seq = self.mutant_seqs[index][0]


		RNA_seq_loop_mut_index = [] # in RNA indexing always
		other_mutations = [] #this is in full indexing if compare_to_full_WT, otherwise RNA indexing
		for mut in self.mutant_seqs[index][1]: #mutant indices
			if self.compare_to_full_WT:
				# convert full index back to RNA index
				mut_ind = mut - (len(self.wt_seq_full) - len(self.wt_RNA_seq))
			else: mut_ind = mut
			if mut_ind > -1:
				# then it's an RNA mutation, compare to the S.S.
				if self.wt_RNA_secstruct[mut_ind] == '.':
					RNA_seq_loop_mut_index.append( mut_ind )
				else:
					# this is a mutation within a S.S. element 
					other_mutations.append( mut )
			else: #this is a protein mutation
				other_mutations.append( mut )
		
		if len( RNA_seq_loop_mut_index ) == 0: # there are no loop mutations to be built w/ stepwise
			# Then just run the low-res protocol
			final_ddg_commands = self.get_low_res_command_lines( index )
			return '', model_cmd_files_with_dir, final_ddg_commands, []
		# Only continuing here if there are loop mutations that need to be built w/ stepwise
		# First need to thread the other_mutations
		# This is part of the setup command
		setup_cmd = ''
		if len(other_mutations) > 0 :
		
			self.write_mutfile( other_mutations, self.mutant_seqs[index][0], 'none', 'mutfile_initial' )

			# Write the setup commands
			# First mutate_and_score_RNP run to make the initial mutations
			setup_cmd += '%smutate_and_score_RNP -ignore_zero_occupancy false ' %(self.rosetta_prefix)
			setup_cmd += '-s %s -mutfile mutfile_initial ' %( self.wt_struct )
			setup_cmd += '-out_prefix initial_muts_ '
			setup_cmd += '%s > make_initial_muts.log\n' %( self.cmd_opts )
			# copy the output structure to "complex_start_struct.pdb"
			setup_cmd += 'cp complex_lowest_scoring_initial_muts_*.pdb complex_start_struct.pdb\n'

		else: # no other mutations that need to be threaded
			setup_cmd += 'cp %s complex_start_struct.pdb\n' %(self.wt_struct)

		# Now we have a starting structure: complex_start_struct.pdb
		# Renumber this pdb for simplicity (renumber back later)
		protein_len = len(self.wt_seq_full) - len(self.wt_RNA_seq)
		rna_len = len(self.wt_RNA_seq)
		full_len = len(self.wt_seq_full)
		setup_cmd += 'renumber_pdb_in_place.py complex_start_struct.pdb A:1-%d B:%d-%d\n' %(protein_len, protein_len+1, full_len)
		setup_cmd += 'renumber_pdb_in_place.py complex_start_struct.pdb A:1-%d B:%d-%d\n' %(protein_len, protein_len+1, full_len)
		# Make a temp renumbered pdb from the wt struct for manipulation here
		os.system('cp %s temp_wt_complex_start_struct.pdb' %(self.wt_struct))
		os.system('renumber_pdb_in_place.py temp_wt_complex_start_struct.pdb A:1-%d B:%d-%d\n' %(protein_len, protein_len+1, full_len))
		os.system('renumber_pdb_in_place.py temp_wt_complex_start_struct.pdb A:1-%d B:%d-%d\n' %(protein_len, protein_len+1, full_len))

		# Loop through each of the RNA_seq_loop_mut_index, create a "subset"
		# For each new subset, figure out if it overlaps with a previous subset
		# if it does, then merge the two subsets
		# a subset should be a set of resnum+residue, i.e. { '1g', '2c', '5u' }
		# chain is unneccessary b/c of renumbering above
		radius = 10 # TODO: Make an option? Change this size?
		mut_subsets = [] # Will be a list of lists of sets of residue numbers and rebuild residues, 1-based indexing, full-complex numbering
		for i in RNA_seq_loop_mut_index:
			rebuild_res = ['B'+ str(i+protein_len+1)]
			rebuild_res_nums = [ i+protein_len+1 ]
			# If we want to rebuild the surrounding residues as well ( for now just resample - 
			# can figure that out from the rebuild_res_nums later )
			#if i+protein_len > protein_len:
			#	rebuild_res.append('B'+str(i+protein_len-1))
			#	rebuild_res_nums.append(i+protein_len-1)
			#if i+protein_len < full_len-2:
			#	rebuild_res.append('B'+str(i+protein_len+1))
			#	rebuild_res_nums.append(i+protein_len+1)
			print rebuild_res
			surr_res, surr_chains = get_surrounding_res.get_surrounding_res('temp_wt_complex_start_struct.pdb', rebuild_res, radius)
			subset = [set(surr_res), rebuild_res_nums]
			subset[0].add(rebuild_res_nums[0])
			# compare this subset to the other subsets in mut_subsets
			mut_subsets_copy = deepcopy( mut_subsets )
			for sub in mut_subsets_copy:
				# if there is any overlap, combine them
				if len(subset[0]&sub[0]) > 0:
					for res in sub[0]:
						subset[0].add(res)
					for rebuild_res in sub[1]:
						subset[1].append(rebuild_res)
					mut_subsets.remove(sub)
			mut_subsets.append(subset)

		# Now add residue identities to the subsets
		if self.compare_to_full_WT:
			mutant_seq_full = self.mutant_seqs[index][0]
		else:
			mutant_seq_full = self.wt_seq_full[0:protein_len] + self.mutant_seqs[index][0]
		mut_subsets_with_ID = []
		for sub in mut_subsets:
			new_res_set = set()
			for resnum in sub[0]:
				new_res_set.add(str(resnum)+mutant_seq_full[resnum-1])
			new_rebuild_res_list = []
			for resnum in sub[1]:
				new_rebuild_res_list.append(str(resnum)+mutant_seq_full[resnum-1])
			mut_subsets_with_ID.append([new_res_set, new_rebuild_res_list])

		# Make a deep copy of the mut_subsets_with_ID to iterate over
		mut_subsets_with_ID_copy = deepcopy( mut_subsets_with_ID )
		mut_depends = [] # List of lists: [mutant index, subset index]


		# Check through the master list of subsets (for all mutants) and see if any of these matches
		# something we're trying to build here, if it does, then delete that subset from the current 
		# list (we don't want to build the same thing twice) and add the mutant index and the subset
		# index to the list of dependencies for this mutant
		for curr_sub in mut_subsets_with_ID_copy:
			for sub_list in self.mutant_subsets:
				for sub in range(len(sub_list[0])):
					#compare the residues sets of the 2 subsets
					if len(sub_list[0][sub][0]-curr_sub[0])==0:
						mut_subsets_with_ID.remove(curr_sub)
						mut_depends.append([sub_list[1],sub])

		
		# Now add these subsets to the master list of subsets for each mutant
		self.mutant_subsets.append([mut_subsets_with_ID,index])

		# Now write the setup cmd
		# Make the threaded "native" structure (only need one of these per mutant
		setup_cmd += '%srna_thread -s complex_start_struct.pdb -o threaded_mut.pdb -seq %s\n' %(self.rosetta_prefix,mutant_seq_full)
		# The resulting pdb needs to be renumbered
		setup_cmd += 'renumber_pdb_in_place.py threaded_mut.pdb A:1-%d B:%d-%d\n' %(protein_len, protein_len+1, full_len)
		setup_cmd += 'renumber_pdb_in_place.py threaded_mut.pdb A:1-%d B:%d-%d\n' %(protein_len, protein_len+1, full_len)


		final_cmd = ''
		# Get rid of the residue IDs in mut_subsets_with_ID, definitely a better way to do this
		mut_subsets_min = [] #minimal set of mut_subsets that needs to be built (equal to mut_subsets_with_ID, without ID)
		for sub_list in mut_subsets_with_ID:
			new_mut_set = set()
			new_mut_list = []
			for r in sub_list[0]: # the set of residues in the subset
				new_mut_set.add( int(r[:-1]) )
			for r in sub_list[1]: # the list of residues being mutated
				new_mut_list.append( int(r[:-1]) )
			mut_subsets_min.append([new_mut_set, new_mut_list])
		
		# Make all the subsets and stepwise command lines
		subset_setup_cmds, subset_model_cmd_files, subset_final_cmds = self.make_high_res_subsets_and_cmd_lines( index, mut_subsets_min, mutant_seq_full, 'threaded_mut.pdb', mut_depends )
		setup_cmd += subset_setup_cmds
		final_cmd += subset_final_cmds
		
		for fil in subset_model_cmd_files:
			newfil = '%d/%s' %(index, fil)
			model_cmd_files_with_dir.append(newfil)
		
		################
		## Concatenate into full structure from subsets (include those from dependencies)
		## For now actually do the combinatorial recombination
		## TODO: if for some reason there are actually a lot of subsets, then this is going to
		## explode pretty quickly 
		## TODO: NEED TO USE mut_depends also!!!!
		#final_cmd += 'rm combined_complex_*.pdb\n' # in case some already exist
		#final_cmd += 'rm top_complexes.out\n' # in case the silent already exists (don't want to append here)
		#for subset_perm in permutations(range(1,21),len(mut_subsets_with_ID)+len(mut_depends)):
		#	final_cmd += 'combine_pdbs -f complex_start_struct.pdb -s '
		#	for i in range(len(mut_subsets_min)):
		#		final_cmd += 'subset_%d.out.%d.pdb ' %(i,subset_perm[i])
		#	for i in range(len(mut_subsets_min),len(mut_depends)):
		#		final_cmd += '../%d/subset_%d.out.%d.pdb ' %(mut_depends[i][0], mut_depends[i][1], subset_perm[i])
		#	pdb_tag = ''
		#	for s in subset_perm:
		#		pdb_tag += str(s)
		#	final_cmd += '-o combined_complex_%s.pdb \n' %(pdb_tag)

		## Combine the pdbs back into a silent file
		#final_cmd += 'combine_silent -s combined_complex_*.pdb -out:file:silent_struct_type binary '
		#final_cmd += '-out:file:silent top_complexes.out\n'

		## Make the mutfile
		#with open('mutfile', 'w') as mfil:
		#	mfil.write('total 0\n-0\n')
		#	mfil.write('input_file top_complexes.out\n')

		## Now do the final mutate_and_score_RNP run to get the complex scores
		#final_cmd += 'mutate_and_score_RNP -ignore_zero_occupancy false '
		#final_cmd += '-s %s -mutfile mutfile ' %( self.wt_struct )
		#final_cmd += '-out_prefix ddG_ '
		#final_cmd += '%s > get_complex_scores.log\n' %( self.cmd_opts )
		############


		# also need to make a list of just mutated positions, for making WT structs
		# already have this: RNA_seq_loop_mut_index
		wt_index = -1
		# Figure out if the appropriate WT structure has already been built
		for wt in range(len(self.high_res_WT_rebuild)):
			if set(RNA_seq_loop_mut_index) == set(self.high_res_WT_rebuild[wt]):
				wt_index = wt
				break
		# If it hasn't been built, then build it here
		if wt_index == -1: # hasn't already been built
			mut_dir = os.getcwd()
			os.chdir('..')
			base_dir = os.getcwd()
			if not os.path.exists('wt_%d' %(len(self.high_res_WT_rebuild)) ):
				os.mkdir('wt_%d' %(len(self.high_res_WT_rebuild)) )
			os.chdir( 'wt_%d' %(len(self.high_res_WT_rebuild)) )
			# Copy the wt start struct here as the "native"
			os.system('cp %s complex_start_struct.pdb\n' %(self.wt_struct))
			# Renumber this for simplicity 
			os.system('renumber_pdb_in_place.py complex_start_struct.pdb A:1-%d B:%d-%d' %(protein_len, protein_len+1, full_len))
			# Try twice (in case some issue, like on biox3)
			os.system('renumber_pdb_in_place.py complex_start_struct.pdb A:1-%d B:%d-%d' %(protein_len, protein_len+1, full_len))
			# Make all the subsets and stepwise command lines
			# TODO: Might be worth figuring out WT subset dependencies (for now I don't think it will actually end up making a difference)
			wt_setup_cmds, wt_model_cmd_files, wt_final_cmds = self.make_high_res_subsets_and_cmd_lines( index, mut_subsets, self.wt_seq_full, 'complex_start_struct.pdb', [])
			#wt_setup_cmds, model_cmd_files, wt_final_cmds = self.make_high_res_subsets_and_cmd_lines( mut_subsets_min, self.wt_seq_full, 'complex_start_struct.pdb', mut_depends)
			setup_cmd += 'cd %s/wt_%d\n' %(base_dir, len(self.high_res_WT_rebuild))
			final_cmd += 'cd %s/wt_%d\n' %(base_dir, len(self.high_res_WT_rebuild))
			setup_cmd += wt_setup_cmds
			final_cmd += wt_final_cmds
			for fil in wt_model_cmd_files:
				newfil = 'wt_%d/%s' %(len(self.high_res_WT_rebuild),fil)
				model_cmd_files_with_dir.append(newfil)
			# Still need to do all the finalization (combination, etc.)

			os.chdir( mut_dir )
			wt_index = len(self.high_res_WT_rebuild)
			self.high_res_WT_rebuild.append(RNA_seq_loop_mut_index)
	
		# write a wt_comparison.txt file in the dir that says which wt struct should be
		# used for comparison
		with open( 'wt_comparison.txt', 'w') as cfil:
			cfil.write( 'wt_' + str(wt_index) + '\n' )

		# List of dependencies for this index, just return a list of the other mutant indices
		# (not specific subsets for simplicity)
		mut_depend_indices = []
		for m in mut_depends:
			mut_depend_indices.append(m[0])

		return setup_cmd, model_cmd_files_with_dir, final_cmd, mut_depend_indices

				
	def make_high_res_subsets_and_cmd_lines( self, index, mut_subsets, mutant_seq_full, native_pdb, mut_depends ):
		setup_cmd = ''
		model_cmd_files = []
		final_cmd = ''
		protein_len = len(self.wt_seq_full) - len(self.wt_RNA_seq)
		full_len = len(self.wt_seq_full)
		for sub in range(len(mut_subsets)):
			# The residues that are being rebuilt: mut_subsets[sub][1] (a list)
			# The residues that are part of the subset: mut_subsets[sub][0] (a set)
			subset_residues = deepcopy( mut_subsets[sub][0] )
			for r in mut_subsets[sub][1]:
				subset_residues.remove( r ) #don't want the mutated residue to be in the input pdb
			subset_residues_string = ' '.join( [str(x) for x in sorted(subset_residues)] )
			setup_cmd += 'pdbslice.py complex_start_struct.pdb -subset %s subset_%d_\n' %(subset_residues_string, sub)
			setup_cmd += 'pdbslice.py complex_start_struct.pdb -subset %s subset_%d_\n' %(subset_residues_string, sub)
			# Make the fasta files
			# loop through all the subset residues
			# Make the tag (need a list of chains first)
			chains = []
			fasta_seq = ''
			for r in sorted(mut_subsets[sub][0]):
				if r < protein_len:
					chains.append('A')
				else:
					chains.append('B')
				fasta_seq += mutant_seq_full[r-1]  # 0 based indexing of seq vs 1 based for subset
			fasta_tag = make_tag_with_dashes( sorted(mut_subsets[sub][0]), chains )
			with open( 'subset_%d.fasta' %(sub), 'w' ) as ffil:
				ffil.write( '> subset_%d_complex_start_struct.pdb %s\n' %(sub,fasta_tag))
				ffil.write( fasta_seq + '\n' )
				
			# Make the stepwise command line, write to README_STEPWISE
			with open('README_STEPWISE_sub_%d' %(sub), 'w') as cmdfil:
				cmdfil.write('%sstepwise -fasta subset_%d.fasta -s subset_%d_complex_start_struct.pdb ' %(self.rosetta_prefix,sub,sub))
				# TODO: is this a good number of cycles?
				# TODO: Maybe this should be based on length of sample_res?
				# Number of stepwise cycles will be equal to 1000*(# residues being built from scratch)
				N_rebuild = len(mut_subsets[sub][1])
				cycles = 1000*N_rebuild
				cmdfil.write( '-virtualize_free_moieties_in_native false -nstruct 100 -cycles %d ' %(cycles) )
				# TODO: assumes there is a variant of the sfxn called sfxn-stepwise.wts
				cmdfil.write( '-score:weights %s-stepwise ' %(self.sfxn)) 
				cmdfil.write( '-out:file:silent subset_%d.out -protein_prepack false ' %(sub) )
				cmdfil.write( '-global_optimize false -switch_focus_frequency 0 -skip_preminimize ' )
				cmdfil.write( '-pack_protein_side_chains false -rmsd_screen 6.0 -native %s ' %(native_pdb))
				# get the sample res:
				sample_res_string = ''
				for rebuild in mut_subsets[sub][1]: # the residues that are being rebuilt
					res = 'B:' + str(rebuild) +' '
					if res not in sample_res_string:
						sample_res_string += res
					if rebuild > protein_len+1:
						res = 'B:' + str(rebuild-1) +' '
						if res not in sample_res_string:
							sample_res_string += res
					if rebuild < (full_len):
						res = 'B:' + str(rebuild+1)+' '
						if res not in sample_res_string:
							sample_res_string += res

				cmdfil.write( '-sample_res %s \n' %(sample_res_string) )

			# Write the submission command line

			if self.run_local:
				setup_cmd += 'source README_STEPWISE_sub_%d\n' %(sub)
			#else:
			# TODO: option to change how long this runs?
			#	setup_cmd += 'rosetta_submit.py README_STEPWISE_sub_%d stepwise_out_sub_%d 50 24\n' %(sub, sub)
			#	setup_cmd += 'source ./qsubMINI\n'
			model_cmd_files.append('README_STEPWISE_sub_%d' %(sub))
			
			# Write the final cmd
			# Collect the results from the stepwise runs
			if not self.run_local:
				final_cmd += 'easy_cat.py stepwise_out_sub_%d\n' %(sub)
			# Extract the low scoring pdbs
			final_cmd += 'extract_lowscore_decoys.py subset_%d.out 20\n' %(sub)
			# Renumber the pdbs back to original numbering
#			# This is not the right place to renumber!! (fixed 03-20)
#			for i in range(1,21):
#				final_cmd += 'renumber_pdb_in_place.py subset_%d.out.%d.pdb %s \n' %(sub,i,self.wt_pdb_tag)
#				final_cmd += 'renumber_pdb_in_place.py subset_%d.out.%d.pdb %s \n' %(sub,i,self.wt_pdb_tag)
#
			###################

		###############
		# Concatenate into full structure from subsets (include those from dependencies)
		# For now actually do the combinatorial recombination
		# TODO: if for some reason there are actually a lot of subsets, then this is going to
		# explode pretty quickly 
		# TODO: NEED TO USE mut_depends also!!!!
		final_cmd += 'rm combined_complex_*.pdb\n' # in case some already exist
		final_cmd += 'rm top_complexes.out\n' # in case the silent already exists (don't want to append here)
		for subset_perm in permutations(range(1,21),len(mut_subsets)+len(mut_depends)):
			# TODO: path for combine_pdbs.py!!!
			final_cmd += 'python ~/bin/combine_pdbs.py -f complex_start_struct.pdb -s '
			for i in range(len(mut_subsets)):
				final_cmd += 'subset_%d.out.%d.pdb ' %(i,subset_perm[i])
			for i in range(len(mut_subsets),len(mut_depends)):
				final_cmd += '../%d/subset_%d.out.%d.pdb ' %(mut_depends[i][0], mut_depends[i][1], subset_perm[i])
			pdb_tag = ''
			for s in subset_perm:
				pdb_tag += str(s)
			final_cmd += '-o combined_complex_%s.pdb \n' %(pdb_tag)

			# Renumber the pdbs back to original numbering
			final_cmd += 'renumber_pdb_in_place.py combined_complex_%s.pdb %s\n' %(pdb_tag, self.wt_pdb_tag)
			final_cmd += 'renumber_pdb_in_place.py combined_complex_%s.pdb %s\n' %(pdb_tag, self.wt_pdb_tag)
		

		# Combine the pdbs back into a silent file
		final_cmd += 'combine_silent -s combined_complex_*.pdb -out:file:silent_struct_type binary '
		final_cmd += '-out:file:silent top_complexes.out\n'

		# Make the mutfile
		with open('mutfile', 'w') as mfil:
			mfil.write('total 0\n-0\n')
			mfil.write('input_file top_complexes.out\n')

		# Now do the final mutate_and_score_RNP run to get the complex scores
		final_cmd += '%smutate_and_score_RNP -ignore_zero_occupancy false ' %(self.rosetta_prefix)
		final_cmd += '-s %s -mutfile mutfile ' %( self.wt_struct )
		final_cmd += '-out_prefix ddG_ '
		pack_reps = self.figure_out_protein_pack_reps( index )
		final_cmd += '-protein_pack_reps %d ' %(pack_reps)
		final_cmd += '%s > get_complex_scores.log\n' %( self.cmd_opts )
		###########

		# Write the REBUILD_RESIDUES.txt file so that neighbors can be repacked and min'ed later
		all_rebuild_residues = []
		for subset in mut_subsets:
			for r in subset[1]: 
				# mut_subsets should be in full complex numbering and should have 1-based indexing
				all_rebuild_residues.append(str(r))
		with open('REBUILD_RESIDUES.txt', 'w') as resfile:
			resfile.write(' '.join(all_rebuild_residues)+'\n')

		return setup_cmd, model_cmd_files, final_cmd
	

	def get_command_lines( self, index ):
		# Make a directory and change into that directory
		curr_dir = os.getcwd()
		if not os.path.exists( curr_dir + '/' + str(index) ):
			os.mkdir( curr_dir + '/' + str(index) )
		os.chdir( curr_dir + '/' + str(index) )

		if self.method == 'low-res':
			final_ddg_commands = self.get_low_res_command_lines( index )
			# there will not be any dependencies for low-res method
			depend_indices = []
			# there are no setup_commands for low-res method
			setup_commands = ''
			# there aren't any stepwise or farfar jobs that need to be run
			model_cmd_files = []
		elif self.method == 'med-res':
			setup_commands, model_cmd_files, final_ddg_commands, depend_indices = self.get_med_res_command_lines( index )
		else: #'high-res'
			setup_commands, model_cmd_files, final_ddg_commands, depend_indices = self.get_high_res_command_lines( index )
		

		os.chdir( curr_dir )

		return setup_commands, model_cmd_files, final_ddg_commands, depend_indices


