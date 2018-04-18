import os
import csv
from itertools import islice
import argparse
import numpy as np
from sys import exit
import glob

wt_ss = '(((((.((....)))))))'
wt_seq = 'acaugaggaucacccaugu'
#wt_seq = 'ACAUGAGGAUCACCCAUGU'

def main():

  seq_index = 0
  run_dir = os.getcwd()
  method = ''
  if args.farfar: 
	method = 'farfar'
  elif args.stepwise:
	if args.stepwise_rna_only:
		method = 'stepwise-rna-only'
	else:
		method = 'stepwise'
	with open('residues_within_10A_of_MS2_chainR.txt', 'r') as resfile:
		protein_res_close = resfile.readline().replace('\n','')
  else:
	print 'Please specify a method for modeling non-canonical base pairs:'
	print 'either -f for farfar or -s for stepwise'
	exit( 0 )

  rebuild_tag = ""
  if args.rebuild_surrounding:
	rebuild_tag+="build_surr"
  else:
	rebuild_tag+="no_build_surr"
  # Create file and directory names
  seq_file_name = "seqs_include_nc_%s_%s%s.txt" %(method,rebuild_tag,args.tag)
  mutfile_name = "%s_%s%s.mutfile" %(method,rebuild_tag,args.tag)
  job_submit_file_name = "submit_%s_%s%s.py" %(method,rebuild_tag,args.tag)
  out_dir_name = "%s_%s%s_add_NC_BPs" %(method,rebuild_tag,args.tag)

  # Don't overwrite files/directories (unless specified)
  j = 0
  if os.path.exists( out_dir_name ) and not args.overwrite:
  	while os.path.exists("%s_%d"%(out_dir_name,j)):
		j+=1
	out_dir_name+="_%d" %(j)
 	seq_file_name+=".%d" %(j)
 	mutfile_name+=".%d" %(j)
 	job_submit_file_name+=".%d" %(j)

  with open ( seq_file_name, 'w') as outfile, open ( mutfile_name, 'w') as mutfile, open(job_submit_file_name, 'w') as job_submit_file:
	job_submit_file.write('import os\n\n')
  	outfile.write( 'Sequence            dG(kT)   RNAdG   #NC  seq_index E_min E_ss p_ss\n' )
  	#with open('small_test.csv') as csvfile:
  	with open('all_MS2_data.csv') as csvfile:
  	   reader = csv.DictReader(csvfile)
	   # Read in the sequence
  	   for row in reader:
		os.chdir( run_dir )
		# Check if there are more than 50 clusters
  		if ( int(row['Clusters']) < 50 ):
			continue
		# Check if we just want to use "high affinity" sequences
		if ( args.high_aff ):
			if ( float(row['\x96_G (kT)']) < 13.13 ):
				continue #go on to the next sequence
		seq = row['Sequence'].lower()
		# Parse base pairs out of secondary structure
		base_pairs = find_pairs( wt_ss )
		# Figure out where there are non-canonical base pairs
		# Default definition of non-canonical is not A-U or G-C
		if ( args.control_canonical ):
			# this means that we want to rebuild any base pair that got mutated
			# so basically want to look through residues that should be base
			# paired (according to the WT secstruct) and see if any residues have
			# been changed
			nc_base_pairs = list_mutated_BPs( base_pairs, seq )
		else:
			nc_base_pairs = list_NC_BPs( base_pairs, seq )
		# Only keep sequences with the specified number of non-canonical base pairs, default=1
		if ( seq != wt_seq ):
			# control_canonical basically changes the meaning of max_NC and num_NC
			# because now *any* mutated base pairs are listed as non-canonical base pairs
			# Here we still want to find out if there are any *true* NC base pairs
			if (args.control_canonical):
				true_nc_base_pairs = list_NC_BPs( base_pairs, seq )
			else: 
				true_nc_base_pairs = nc_base_pairs
			# Check if user specified a max number of base pairs
			if ( args.max_NC is not None ):
				# if so, check if this sequence exceeds the max num allowed NC bps
				if ( len(true_nc_base_pairs) > args.max_NC ):
					continue # go on to the next sequence
			# Otherwise check if there are exactly num_NC non-canonical base pairs
			elif ( len(true_nc_base_pairs) != args.num_NC ): #initial set, default just look at 1 non-canonical
				continue # go on to the next sequence

		if ( args.no_extra_mutations ):
			############################
			# For initial testing, just look at sequences with one non-canonical
			# count mutations
			mut_count = 0
			# Fixing this logic 1/8/16
			nc_residues = []
			for ind in nc_base_pairs:
				nc_residues.append(base_pairs[ind]['start'])
				nc_residues.append(base_pairs[ind]['end'])
			##
			for i in range(len(seq)):
				if i in nc_residues: continue
				#if i in nc_base_pairs: continue # changed on 1/8/16
				if seq[i]!=wt_seq[i]: mut_count+=1; break
			#if ( (seq[5] != 'a') or (seq[8] != 'a') or (seq[9] != 'u') or (seq[10]!='c') or (seq[11]!='a') ):
			if (mut_count>0):
				continue #go on to the next sequence
			#############################

		# Finally, check if we only want to keep sequences that preserve the WT sec struct
		if ( args.maintains_WT_SS ):
			# Then check the secondary structure
			# TODO: this runs RNAfold, but it will be run below too, wasteful!
			pred_ss = get_ss( seq )
			if pred_ss != wt_ss:
				continue # go on to the next sequence

		# Check if there are any non-canonical base pairs
		if ( len(nc_base_pairs) > 0 ):
			# Make a "canonical mutfile"
			# Find all sequence mutations that are not part of a new non-canonical bp
			# "Flatten" dict to list for later simplicity
			list_of_nc_pos = []
			for i in nc_base_pairs:
				list_of_nc_pos.append(base_pairs[i]['start'])
				list_of_nc_pos.append(base_pairs[i]['end'])
			list_surrounding = []
			# Make a new directory for this seq_index
			if not os.path.exists( out_dir_name ):
				os.mkdir( out_dir_name )
			if not os.path.exists( '%s/%d' %(out_dir_name,seq_index) ):
				os.mkdir( '%s/%d' %(out_dir_name,seq_index) )
			os.chdir( '%s/%d' %(out_dir_name,seq_index) )
			with open( 'mutfile_c_%d' %(seq_index), 'w' ) as cmfile:
				mut_string = ''
				total_muts = 0
				if args.rebuild_surrounding:
					# Add surrounding base pairs
					for bp in nc_base_pairs:
						if bp > 0:
							# the base pair below
							list_surrounding.append(base_pairs[bp-1]['start'])
							list_surrounding.append(base_pairs[bp-1]['end'])
						if bp < len(base_pairs)-1:
							# the base pair above
							list_surrounding.append(base_pairs[bp+1]['start'])
							list_surrounding.append(base_pairs[bp+1]['end'])
				for i in range(0,len(seq)):
					# check that it's not a non-canonical (or surrounding?)
					if i in list_of_nc_pos: continue
					#if (i in list_of_nc_pos) or (i in list_surrounding): continue
					if (seq[i] != wt_seq[i]):
						mut_string+=(wt_seq[i] + ' ' + str(i) + ' ' + seq[i] + '\n')
						total_muts+=1
				cmfile.write( 'total ' + str(total_muts) + '\n')
				cmfile.write(str(total_muts) + '\n')
				cmfile.write( mut_string )
			if ( total_muts > 0 ): # only need to run this if there are "canonical mutations" also
				# Now run mutate_and_score app to get new "bound RNA alone" conf
				if (args.start_struct != ""):
					os.system('~/src/Rosetta_second/main_RNP/source/bin/mutate_and_score_RNP -ignore_zero_occupancy false -s ../../%s -score:weights ../../P_overlap_reR-hbond_sc.wts -move_backbone false -out_prefix initial_muts_ -mutfile mutfile_c_%s -unbound_protein ../../test_cutoff200_wildtype_unbound_protein.pdb  -mutate_unbound_RNA true -prot_offset_num 258 -RNA_offset_num 1 -protein_start_res 1 -protein_end_res 258 -relax_cutoff_dist 20.0 -dump_bound_rna true -bound_rna_dump_tag seq_%s >initial_muts.log' %(args.start_struct,str(seq_index),str(seq_index)) )
				else:
					os.system('~/src/Rosetta_second/main_RNP/source/bin/mutate_and_score_RNP -ignore_zero_occupancy false -s ../../min_again_relax_full_RNA_P_overlap_reR-hbond_sc_wildtype_bound.pdb -score:weights ../../P_overlap_reR-hbond_sc.wts -move_backbone false -out_prefix initial_muts_ -mutfile mutfile_c_%s -unbound_protein ../../test_cutoff200_wildtype_unbound_protein.pdb  -mutate_unbound_RNA true -prot_offset_num 258 -RNA_offset_num 1 -protein_start_res 1 -protein_end_res 258 -relax_cutoff_dist 20.0 -dump_bound_rna true -bound_rna_dump_tag seq_%s >initial_muts.log' %(str(seq_index),str(seq_index)) )
				# And wait here until it finishes running, this step is slow...
				#TODO: generalize chain!!!
				if args.farfar:
					# For farfar, get the RNA structure for input into fragment assembly
					initial_muts_start_struct = glob.glob('initial_muts*_bound.pdb')
					os.system('extract_chain.py %s R' %(initial_muts_start_struct[0]))
					os.system('mv %sR.pdb RNA_bound_conf_seq_%d.pdb' %(initial_muts_start_stuct[0].replace('.pdb',''), seq_index))
				# If using stepwise, we want the full bound structure
				if args.stepwise:
					os.system('mv initial_muts_*_bound.pdb stepwise_input_seq_%d.pdb' %(seq_index))
					#TODO: be careful about recombining this later with the full protein, if not included in stepwise run!!!
			else: 
				# copy the wt rna bound structure
				if args.farfar:
					if (args.start_struct != ""):
						# TODO: this is super hacky, generalize this!!
						os.system('extract_chain.py ../../%s R' %(args.start_struct))
						rna_struct = (args.start_struct).replace('.pdb','') + 'R.pdb'
						print 'cp ../../%s RNA_bound_conf_seq_%d.pdb' %(rna_struct, seq_index)
						os.system('pwd')
						os.system('cp ../../%s RNA_bound_conf_seq_%d.pdb' %(rna_struct, seq_index) )
					else:
						os.system('cp ../../min_again_relax_full_RNA_P_overlap_reR-hbond_sc_wildtype_bound_RNA_only.pdb RNA_bound_conf_seq_%d.pdb' %(seq_index) )
				if args.stepwise:
					if (args.start_struct != ""):
						os.system('cp ../../%s stepwise_input_seq_%d.pdb' %(args.start_struct, seq_index) )
					else:
						os.system('cp ../../min_again_relax_full_RNA_P_overlap_reR-hbond_sc_wildtype_bound.pdb stepwise_input_seq_%d.pdb' %(seq_index) )
			if args.farfar:
				#TODO: is this actually necessary?
				# First do a quick renumber to not mess things up 
				os.system('renumber_pdb_in_place.py RNA_bound_conf_seq_%s.pdb R:1-%d' %(str(seq_index),len(seq)))

			# Keep only the protein residues within 10A of RNA, I used pymol to generate the text file I read in here
			# TODO: more general way to get close protein residues
			if args.stepwise:
				if args.stepwise_rna_only: 
					os.system('pdbslice.py stepwise_input_seq_%d.pdb -subset R:1-%d subset_' %(seq_index, len(seq)))
				else:
					os.system('pdbslice.py stepwise_input_seq_%d.pdb -subset %s R:1-%d subset_' %(seq_index,protein_res_close,len(seq)))
			# Now take the new RNA structure and excise the non-canonical positions
			# Figure out which positions to excise
			excise_str = ''
			for i in list_of_nc_pos:
				excise_str+=('R'+str(i+1))
				excise_str+=' '
			# TODO: chain for excise should be general

			##############
			# if any of the last three base pairs is being modified, rebuild it and all the base pairs below it
			# TODO: this needs to be seriously generalized
			end_bp_res = [2,1] # don't need to include first basepair b/c nothing is "below" it
			cut_below_point = 0
			# TODO: doing this seems to cause some stepwise bug!? comment out for now
			for i in end_bp_res: # checking for nc muts in last 3 basepairs
				if base_pairs[i]['start'] in list_of_nc_pos:
					#this means that we're remodeling one of the last 3 non-canonical base pairs
					cut_below_point = i
					break # checking in highest to lowest order, so this is the top of where to cut
			if cut_below_point > 0: 
				for i in range( cut_below_point ):
					excise_str+=('R'+str(base_pairs[i]['start']+1)+' ')
					excise_str+=('R'+str(base_pairs[i]['end']+1)+' ')
			#################
			if args.rebuild_surrounding:
				# If specified, also excise the surrounding base pairs
				for i in list_surrounding:
					excise_str+=('R'+str(i+1)+' ')
			#################

			# Do the excision!
			if args.farfar:
				os.system('pdbslice.py RNA_bound_conf_seq_%s.pdb -excise %s scaff_' %(str(seq_index),excise_str))
				# Sec struct for rna_denovo, can't have non-canonical "base pairs"
				seq_ss_list = list(wt_ss)
				for i in list_of_nc_pos:
					seq_ss_list[i] = '.'
				new_ss_str = ''.join(seq_ss_list)
				with open('seq_%d.fasta' %(seq_index),'w') as fasta:
					fasta.write( '> RNA_bound_conf_seq_%d.pdb R:1-%d\n' %(seq_index,len(seq)))
					fasta.write( '%s' %(seq))
				# Make rna_denovo input files
				os.system('rna_denovo_setup.py -nstruct 100 -fasta seq_%d.fasta -secstruct "%s" -fixed_stems -chemical:enlarge_H_lj -s scaff_RNA_bound_conf_seq_%d.pdb -tag seq_%d -score:weights ../../P_overlap_reR-hbond_sc.wts' %(seq_index, new_ss_str, seq_index, seq_index))
				# This creates a README_FARFAR
				# HACK: change the rosetta directory TODO: fix for real
				os.system('sed -i "s/main/main_RNP/g" README_FARFAR')
				os.system('rosetta_submit.py README_FARFAR farfar_out 1 4')
			
			if args.stepwise:
				os.system('pdbslice.py subset_stepwise_input_seq_%d.pdb -excise %s scaff_' %(seq_index,excise_str))
				# use pdb2fasta.py to get initial fasta which we can modify with correct seq
				os.system('pdb2fasta.py scaff_subset_stepwise_input_seq_%d.pdb > temp.fasta' %(seq_index))
				with open('temp.fasta','r') as fasta, open('seq_%d.fasta' %(seq_index), 'w') as new_fasta:
					for line in fasta:
						# look for the rna sequence line
						if line.startswith('>'):
							new_fasta.write(line)
						elif ('a' in line) or ('u' in line) or ('g' in line) or ('c' in line):
							new_fasta.write('%s' %(seq))
						else: 
							new_fasta.write(line)
				os.remove('temp.fasta')
				# Make stepwise input file, only one input structure
				with open('README_STEPWISE', 'w') as readme:
					readme.write('~/src/Rosetta_second/main_RNP/source/bin/stepwise -fasta seq_%d.fasta -s scaff_subset_stepwise_input_seq_%d.pdb -virtualize_free_moieties_in_native false -skip_preminimize -nstruct 100 -cycles 2000 -score:weights ../../P_overlap_reR-hbond_sc_stepwise.wts -out:file:silent seq_%d.out\n' %(seq_index,seq_index,seq_index))

				os.system('rosetta_submit.py README_STEPWISE stepwise_out 6 48')
				#os.system('rosetta_submit.py README_STEPWISE stepwise_out 1 10')
			
			## Print out a file that lists all the excised residues so that we can read it later in the main app
			with open('REBUILD_RESIDUES.txt','w') as resfile:
				resfile.write(excise_str.replace('R','')+'\n')
			# Add some submission lines to farfar_submit file
			job_submit_file.write('os.chdir("' + run_dir + '")\n')
			job_submit_file.write('os.chdir("%s/' %(out_dir_name) + str(seq_index) + '")\n')
			job_submit_file.write('os.system("source ./qsubMINI")\n')
			# Running this file after running this script will submit everything
			# Add this sequence to the full mutfile
			# There shouldn't actually be any mutations done in the full run
			#total_muts, mut_string = get_mut_string( seq )
			mutfile.write('total 0\n')
			mutfile.write('0\n')
			# This will be the lowest score structure from the farfar run
			mutfile.write('input_file %s/%d/top20_seq_%d.out \n' %(out_dir_name,seq_index,seq_index) )
			#mutfile.write('input_file farfar_add_NC_BPs/%d/seq_%d.out.1.pdb \n' %(seq_index,seq_index) )
		else: # no non-canonical base pairs
			# Add the sequence to the full mutfile
			total_muts, mut_string = get_mut_string( seq )
			mutfile.write('total ' + str(total_muts) + '\n')
			mutfile.write(str(total_muts) + '\n')
			mutfile.write( mut_string )
			# There is no input structure to put in the mutfile
			mutfile.write('input_file none\n')
		# TODO:
		# Write out a qsub file, replace mutfile, starting structure, dir?

		# Figure out the free energy of the sequence (want to write this to outfile)
		ftemp = open('temp_seq.txt', 'w')
		ftemp.write( seq + '\n')
		ftemp.close()
		# Use RNAfold to get the free energy
		p_min = ''
		E_min = ''
		fe = ''
		os.system('RNAfold -p -T 25 <temp_seq.txt> temp_ss.out')
		with open('temp_ss.out', 'r') as ssfile:
			line_index = 0
			for l in ssfile:
				if line_index == 2:
					if len(l.split())>2:
						fe_ = l.split()[2]
						fe = fe_.replace(']','')
					else: 
						fe__ = l.split()[1]
						fe_ = fe__.replace('[','')
						fe = fe_.replace(']','')
				if line_index == 4:
					p_min_ = l.split()[6]
					p_min = p_min_.replace(';','')
				# Get the minimum energy
				if line_index == 1:
					if len(l.split())>2:
						E_min_ = l.split()[2]
						E_min = E_min_.replace(')','')
					else:
						E_min__ = l.split()[1]
						E_min_ = E_min__.replace('(','')
						E_min = E_min_.replace(')','')
				line_index +=1

		# Want to get the probability of wt ss formation as well
		# Get the wildtype secondary structure, with any non-canonical
		# (excluding GU pairs) pairs listed as '.' instead of '('
		nc_base_pairs_noGU = list_NC_BPs_noGU( base_pairs, seq )
		wt_ss_list = list(wt_ss)
		for k in nc_base_pairs_noGU:
			for j in base_pairs[k].values():
				wt_ss_list[j] = '.'
		ss_str_noNC = ''.join(wt_ss_list)
		E_ss = ''
		ftemp = open('temp_seq.txt', 'w')
		ftemp.write( seq + '\n' )
		ftemp.write( ss_str_noNC + '\n' )
		ftemp.close()
		os.system('RNAeval -T 25 <temp_seq.txt> temp_eval.out')
		with open('temp_eval.out', 'r') as evalfile:
			for l in islice(evalfile, 1,2):
				if len(l.split())>2:
					E_ss_ = l.split()[2]
					E_ss = E_ss_.replace(')','')
				else:
					E_ss__ = l.split()[1]
					E_ss_ = E_ss__.replace('(','')
					E_ss = E_ss_.replace(')','')


		# RNAeval with wt ss
		ftemp = open('temp_seq.txt', 'w')
		ftemp.write( seq + '\n' )
		ftemp.write( wt_ss + '\n' )
		ftemp.close()
		os.system( 'RNAeval -T 25 <temp_seq.txt> temp_eval.out')
		with open('temp_eval.out', 'r') as evalfile:
			for l in islice(evalfile,1,2):
				if len(l.split())>2:
					E_wtss_ = l.split()[2]
					E_wtss = E_wtss_.replace(')','')
				else:
					E_wtss__ = l.split()[1]
					E_wtss_ = E_wtss__.replace('(','')
					E_wtss = E_wtss_.replace(')','')

		# Calculate the probability of this wt ss formation:
		kT = 0.59
		p_ss = float(p_min) * np.exp((float(E_min)-float(E_ss))/kT)

		# Write data to the output file
		outfile.write( seq + ' ' + row['\x96_G (kT)'] + ' '+ fe + ' '+ str(len(nc_base_pairs)) + ' ' + str(seq_index) )
		outfile.write ( ' ' + str(E_min) + ' ' + str(E_ss) + ' ' + str(p_ss) + ' ' + E_wtss + '\n' )
		seq_index += 1

  print("Writing %s files out to %s" %(method, out_dir_name))
  print("Submission file is %s" %(job_submit_file_name)) 
  print("Run 'python %s' to submit %s jobs" %(job_submit_file_name, method))
  print("Sequences and nearest neighbor free energies written to %s" %(seq_file_name))
  


####### Some helper functions #######
def get_mut_string( seq ):
	mut_string = ''
	total_muts = 0
	for i in range(0,len(seq)):
		if (seq[i] != wt_seq[i]):
			mut_string+=(wt_seq[i] + ' ' + str(i) + ' ' + seq[i] + '\n')
			total_muts+=1
	return total_muts, mut_string

def check_nostruct_matches( sequence ):
	nostruct_residues = [ 0, 1, 2, 16, 17, 18 ]
	for i in nostruct_residues:
		if (sequence[i] != wt_seq[i]):
			return False
	return True

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

def list_NC_BPs_noGU( bps, seq ):
	list_nc_bps = []
	# loop through the list of base pairs
	for i in range(len(bps)):
		if ((seq[bps[i]['start']] == 'a') and (seq[bps[i]['end']] == 'u')
			or (seq[bps[i]['start']] == 'u') and (seq[bps[i]['end']] == 'a')
			or (seq[bps[i]['start']] == 'g') and (seq[bps[i]['end']] == 'c')
			or (seq[bps[i]['start']] == 'c') and (seq[bps[i]['end']] == 'g')
			or (seq[bps[i]['start']] == 'g') and (seq[bps[i]['end']] == 'u')
			or (seq[bps[i]['start']] == 'u') and (seq[bps[i]['end']] == 'g')):
			continue
		list_nc_bps.append(i)

	# Return a list of indices of nc base pairs, can be used to index bps
	return list_nc_bps 

def find_offset( ss ):
	offsets = []
	offset_list = []
	# brute force it and loop through all possible offsets
	for offset in xrange(-len(wt_ss)+1,len(wt_ss)):
		num_matches = 0
		for i in range(len(wt_ss)):
			if i-offset < len(ss):
				if wt_ss[i] == ss[i-offset]: num_matches+=1
		offsets.append(num_matches)
		offset_list.append(offset)

	# This is the best offset to use
	# if this is 0, then best not to apply an offset
	return offset_list[np.argmax(np.array(offsets))]

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

# Takes a sequence and returns the ss
def get_ss( seq ):
	# make a temporary sequence file
	ftemp = open('tmp_seq.txt', 'w')
	ftemp.write( seq +'\n')
	ftemp.close()
	os.system('RNAfold -T 21 <tmp_seq.txt> tmp_ss.out')
	ss = ''
	with open( 'tmp_ss.out', 'r') as ssfile:
		for line in islice( ssfile, 1, 2): 
			ss = line.split()[0]
	
	# remove the temporary sequence and output files
	os.system('rm tmp_seq.txt')
	os.system('rm tmp_ss.out')
	
	return ss


if __name__ == '__main__':
	print "\n#################################################"
	print "WARNING: Please use updated versions of these scripts in apps/public/rnp_ddg/"
	print "#################################################\n"
	parser=argparse.ArgumentParser(description="Set up ddg calculations for list of RNP mutants")
	parser.add_argument('-f', '--farfar', action='store_true', default=False, help='Use farfar to remodel non-canonical base pairs')
	parser.add_argument('-s', '--stepwise', action='store_true', default=False, help='Use stepwise to remodel non-canonical base pairs')
	parser.add_argument('--stepwise_rna_only', action='store_true', default=False, help="Don't include protein in stepwise remodeling")
	parser.add_argument('-rebuild_surrounding', '--rebuild_surrounding', dest='rebuild_surrounding', action='store_true', help='Rebuild base pairs above and below the non-canonical base pair')
	parser.add_argument('-no_rebuild_surrounding', '--no_rebuild_surrounding', dest='rebuild_surrounding', action='store_false', help='Rebuild base pairs above and below the non-canonical base pair')
	parser.set_defaults(rebuild_surrounding=False)
	parser.add_argument('-t', '--tag', type=str, default="", help='Tag for naming files and directories')
	parser.add_argument('-w', '--overwrite', action='store_true', default=False, help='Overwrite directory with default name if it exists')
	parser.add_argument('--start_struct', type=str, default="", help='Starting complex structure to use')
	parser.add_argument('--high_aff', action='store_true', default=False, help='Only look at sequences with affinity below -13.13kT')
	parser.add_argument('--num_NC', type=int, default=1, help='Keep sequences with only this number of non-canonical bps')
	parser.add_argument('--control_canonical', action='store_true', default=False, help='Rebuild canonical base pair mutations with farfar or stepwise as a control')
	parser.add_argument('--no_extra_mutations', action='store_true', default=False, help='Allow only non-canonical base pair mutations')
	parser.add_argument('--maintains_WT_SS', action='store_true', default=False, help='Keep only seqs predicted to have WT sec struct')
	parser.add_argument('--max_NC', type=int, help='Keep only sequences with fewer than this number of non-canonical base pairs')
	args = parser.parse_args()
	main()
