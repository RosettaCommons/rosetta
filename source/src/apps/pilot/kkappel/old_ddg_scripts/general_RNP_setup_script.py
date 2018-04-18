import os
import argparse
import numpy as np
import subprocess
import mutant_modeler


# Do output setup (make informative tag from options, don't overwrite directories)
# make a setup directory, each sequence should get its own directory within
# this setup directory
# - will need to figure out from the input struct which chains are RNA and which are protein
# To be run in the setup directory:
# - do all threading first, if only threading needs to be done (low-res), can still just do it in this dir
# - next do either stepwise loop remodeling (high-res) or farfar bps (medium-res)
# - thread canonical base pair mutations, if there are other non-canonical base pair mutations or loop residue mutations
# - non-canonical base pair mutation remodeling (for now only do with farfar or threading)
# - loop residues mutation remodeling
# - how to get back to the full structure...what other mutations need to be included

# "Finalize" step, what will need to be done: (this will all be written into a qsub file)
# - stepwise structures will need to be reassembled into full structs from subsets
# - farfar threaded structures will need to be reassembled with protein components
# - each directory should end up with a silent file (or a link to a silent file - for the low-res stuff should be a path to silent file containing many relaxed structs (?))
# - a mutfile (can be written initially)
# - a qsub file for the ddG calculation (can be written initially) -- this qsub file should just contain the setup/cleanup commands prior to the mutate_and_score command 
# - if I do this really intelligently, this last qsub file could be submitted by this script with a job dependency (might be a little tricky)
# qsub job can have multiple dependencies: -W depend=afterok:12345:12346 (just colon separate the nums)

# Then need one final script that parses ddG output and gets complex scores, reads in nearest neighbor energies and then spits out a file with predicted ddGs

# if object oriented, might want "LowResModeler", "MedResModeler", "HighResModeler"
# Modeler
# get_reusable_subset
# stores all subsets that have already been built
# need data structure like, subset (moveable residues, fixed residues), jobID, seq#

# all mutate_and_score_RNP options should be passed in to this script
def main( args ):
	# Figure out what method was specified, default: low-res (threading only)
	if args.high_res: method = 'high-res'
	elif args.med_res: method = 'med-res'
	else: method = 'low-res'

	# setup the output file and directory names
	# write everything in this directory, then all files can have standard names
	output_dir = setup_file_and_dir_names( args.tag, args.overwrite, method )

	# Get the absolute path of the start_struct to avoid trouble later
	start_struct_path = os.path.abspath( args.start_struct )

	# Set up the mutant modeler, this will get the wt sequence (full and RNA)
	# and it will figure out the RNA and protein chains
	modeler = mutant_modeler.MutantModeler( method, start_struct_path )
	if args.wt_secstruct != '': 
		modeler.set_wt_RNA_secstruct( args.wt_secstruct )
	if args.submit_qsub_files:
		modeler.set_run_locally(False)
	modeler.set_sfxn( args.sfxn )
	modeler.set_rosetta_prefix( args.rosetta_prefix )

	# If there are additional command line options, add them to the modeler
	cmd_opts = '-score:weights %s ' %( args.sfxn )
	cmd_opts += '-relax_cutoff_dist %s ' %( args.relax_cutoff_dist )
	modeler.set_protein_pack_reps( args.protein_pack_reps )
	cmd_opts += '-Nreps %s ' %( args.Nreps )
	if args.no_min_jumps: min_jumps = 'false'
	else: min_jumps = 'true'
	if args.move_backbone: move_bb = 'true'
	else: move_bb = 'false'
	cmd_opts += '-min_jumps %s ' %( min_jumps )
	cmd_opts += '-move_backbone %s ' %( move_bb )
	if args.move_protein_backbone:
		cmd_opts += '-move_protein_backbone true '
	# TODO: HACK FOR TESTING ONLY
	#cmd_opts += '-min_only true '
	
	modeler.set_cmd_opts( cmd_opts )

	# Go through all the mutants and add them to the modeler
	mutant_indices = [0] # 0 is the WT seq
	with open( args.seq_file, 'r' ) as sfil:
		for line in sfil:
			mut_index = modeler.add_mutant( line.replace('\n','') )
			mutant_indices.append( mut_index )
	
	# Make the base directory and change into it so all new files will be written there
	script_run_dir = os.getcwd()
	if not os.path.exists( script_run_dir + '/' + output_dir ):
		os.mkdir( script_run_dir + '/' + output_dir )
	os.chdir( script_run_dir + '/' + output_dir )

	# Write out the settings that were used to run this script
	with open('general_setup_settings.txt', 'w') as sfil:
		sfil.write('METHOD: %s\n' %(method))
		sfil.write('TAG: %s\n' %(args.tag))
		sfil.write('OVERWRITE?: %s\n' %(str(args.overwrite)))
		sfil.write('START_STRUCT: %s\n' %(start_struct_path))
		sfil.write('WT_SECSTRUCT (user specified): %s\n' %(args.wt_secstruct))
		sfil.write('WT_SECSTRUCT: %s\n' %(modeler.wt_RNA_secstruct))
		sfil.write('SEQ_FILE: %s\n' %(args.seq_file))
		sfil.write('Extra command line options: %s\n' %(cmd_opts))
		

	# TODO: Might want to reorder list of mutants by number of mutations, in some cases might save computation

	local_file = open('all_commands.txt', 'w')
	job_IDs = {}
	# Get command lines for all the mutants
	for mutant_index in mutant_indices:
		setup_commands, model_cmd_files, final_ddg_commands, depend_indices = modeler.get_command_lines( mutant_index )
		# Make qsub files, if specified	
		run_dir = script_run_dir + '/' + output_dir + '/' + str(mutant_index)
		if args.submit_qsub_files: 
			# Only need to write a qsub_initial file if there is a setup command
			setup_jid = ''
			if len(setup_commands) > 0:
				with open('%d/qsub_initial' %(mutant_index), 'w') as qifil:
					# Write qsub header
					qifil.write( '#!/bin/bash\n')
					qifil.write('#PBS -o %s/job_initial.out\n' %(run_dir))
					qifil.write('#PBS -e %s/job_initial.err\n' %(run_dir))
					qifil.write('#PBS -m n\n')
					qifil.write('#PBS -M nobody@stanford.edu\n')
					qifil.write('#PBS -l walltime=24:00:00\n')
					qifil.write('#PBS -N setup_mut_%d\n\n' %(mutant_index))
					# No dependencies for setup commands
					qifil.write('cd %s \n\n' %(run_dir))
					qifil.write(setup_commands + '\n') # For med-res and high-res, this does not include the last qsubMINI submission
				# Submit the setup job and get the jid
				try:
					output = subprocess.check_output('qsub %d/qsub_initial' %(mutant_index), shell=True)
				except: 
					output = subprocess.check_output('qsub %d/qsub_initial' %(mutant_index), shell=True)
				setup_jid = output.split('.')[0]
				print setup_jid
			# Setup the jobs that will run farfar or stepwise, from the model_cmd_files
			remodel_jid = ''
			if method=='med-res' and len(model_cmd_files) > 0:
				# Write a qsub file that runs source README_FARFAR for each README_FARFAR (may be mutant and WT)
				for readme in model_cmd_files:
					fil_dir = readme.split('/')[0]
					os.chdir(fil_dir)
					run_d = os.getcwd()
					with open('qsub_remodel', 'w') as qrm:
						qrm.write('#!/bin/bash\n')
						qrm.write('#PBS -o %s/job_remodel.out\n' %(run_d))
						qrm.write('#PBS -e %s/job_remodel.err\n' %(run_d))
						qrm.write('#PBS -m n\n')
						qrm.write('#PBS -M nobody@stanford.edu\n')
						qrm.write('#PBS -N remodel_mut_%s\n' %(fil_dir))
						qrm.write('#PBS -l walltime=24:00:00\n') # 24 hours should be enough for any number of structs/cycles
						# Depends on the setup job
						qrm.write('#PBS -W depend=afterany:%s\n' %(setup_jid))
						qrm.write('cd %s \n\n' %(run_d))
						qrm.write('source ./%s\n' %(readme.split('/')[-1])) #TODO: might want to actually get this file from model_cmd_files
					# Submit the remodel job and get the jid
					try:
						output = subprocess.check_output('qsub qsub_remodel', shell=True)
					except:
						output = subprocess.check_output('qsub qsub_remodel', shell=True)
					#remodel_jid = output.split('.')[0]
					if remodel_jid == '':
						remodel_jid += '%s' %(output.split('.')[0])
					else:
						remodel_jid += ':%s' %(output.split('.')[0])
					os.chdir('..')
			if method=='high-res' and len(model_cmd_files) > 0: 
				# for each README_STEPWISE_* file, run rosetta_submit.py
				for readme in model_cmd_files:
					print "working with %s" %(readme)
					# open the file and figure out how many sample res there are to determine how many jobs to submit
					with open(readme, 'r') as rfil:
						print "opened the file"
						for opts in rfil.readlines():
							for opt in opts.split('-'):
								if "sample_res" in opt:
									sample_res = opt.split(' ')
									while '' in sample_res: sample_res.remove('')
									while '\n' in sample_res: sample_res.remove('\n')
									nres = len(sample_res) - 1
									break
					# Submit jobs
					# go into the correct directory
					fil_dir = readme.split('/')[0]
					os.chdir(fil_dir) #either mutant index or some wt_ dir
					njobs = nres*10
					out_suffix = readme.split('_')[-1]
					#out_suffix = readme.replace('README_STEPWISE_', '')
					os.system('rosetta_submit.py %s stepwise_out_sub_%s %d 48' %(readme.split('/')[-1], out_suffix, njobs))
					# Add dependencies to the qsubMINI submission file
					os.system('sed -i "s/qsub /qsub -Wdepend=afterany:%s /g" qsubMINI' %(setup_jid))
					#os.system('sed -i "s/qsub /qsub -Wdepend=afterok:%s /g" qsubMINI' %(setup_jid))
					# Submit the jobs
					# Assuming jobs will run sequentially, so only use last jid as remodel_jid for next dependency
					# TODO: this is a very stupid assumption...
					output = subprocess.check_output('source ./qsubMINI', shell=True)
					jid = output.split('\n')[len(output.split('\n'))-2].split('.')[0]
					if remodel_jid == '':
						remodel_jid += '%s' %(jid)
					else:
						remodel_jid += ':%s' %(jid)
					os.chdir('..')
			# Write the final qsub file, cleans up from setup and remodel jobs, and does the final mutate_and_score run to get ddG		
			with open('%d/qsub_final' %(mutant_index), 'w') as qffil:
				# Write qsub header
				qffil.write( '#!/bin/bash\n')
				qffil.write('#PBS -o %s/job_final.out\n' %(run_dir))
				qffil.write('#PBS -e %s/job_final.err\n' %(run_dir))
				qffil.write('#PBS -m n\n')
				qffil.write('#PBS -M nobody@stanford.edu\n')
				qffil.write('#PBS -l walltime=24:00:00\n')
				qffil.write('#PBS -N final_mut_%d\n' %(mutant_index))

				# Add dependencies 
				qffil.write('#PBS -W depend=afterany') #TODO: might actually want afterany
				if setup_jid != '': 
					qffil.write(':%s' %(setup_jid))
				if remodel_jid != '':
					qffil.write(':%s' %(remodel_jid))
				for depend in depend_indices:
					qffil.write(':%s' %(job_IDs[depend]))
				qffil.write('\n\n')
				
				qffil.write('cd %s \n\n' %(run_dir))
				qffil.write(final_ddg_commands)
			# Then submit the qsub files
			try: 
				output = subprocess.check_output('qsub %d/qsub_final' %(mutant_index), shell=True)
			except:
				output = subprocess.check_output('qsub %d/qsub_final' %(mutant_index), shell=True)
			final_jid = output.split('.')[0]
			job_IDs[mutant_index]=final_jid
			
		# Write to a single file that can be run locally
		# Change to the appropriate directory
		local_file.write('cd ' + script_run_dir + '/' + output_dir + '/' + str(mutant_index) + '\n')
		local_file.write(setup_commands + '\n')
		local_file.write('cd ' + script_run_dir + '/' + output_dir + '/' + str(mutant_index) + '\n')
		local_file.write(final_ddg_commands + '\n')
		# Write nearest neighbor energy


		# Figure out the free energy of the sequence (want to write this to outfile)
		ftemp = open('temp_seq.txt', 'w')
		ftemp.write( modeler.mutant_seqs[mutant_index][0] + '\n')
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
		os.remove('temp_seq.txt')
		with open('%d/rna_nn_E.txt' %(mutant_index), 'w') as nnfil:
			nnfil.write('FE   E_min   p_min\n')
			nnfil.write('%s  %s  %s\n' %(fe, E_min, p_min))


	local_file.close()

def setup_file_and_dir_names( tag, overwrite, method ):
	base_out_dir = 'ddG_' + tag + '_' + method
	# Don't overwrite files/directories (unless specified)
	if os.path.exists( base_out_dir ) and not overwrite:
		j = 0
		while os.path.exists("%s_%d" %(base_out_dir, j)):
			j+=1
		base_out_dir+="_%d" %(j)
	return base_out_dir

if __name__ == '__main__':
	print "\n#################################################"
	print "WARNING: Please use updated versions of these scripts in apps/public/rnp_ddg/"
	print "#################################################\n"
	parser=argparse.ArgumentParser(description="Set up ddg calculations for list of RNP mutants")
	parser.add_argument('--high_res', action='store_true', default=False, help='Make loop mutations using stepwise, all others threaded')
	parser.add_argument('--med_res', action='store_true', default=False, help='Make non-canonical base pair mutations using farfar, all others threaded')
	parser.add_argument('--low_res', action='store_true', default=False, help='Make all mutations with simple threading protocol (DEFAULT METHOD)')
	parser.add_argument('-t', '--tag', type=str, default="", help='Tag for naming files and directories')
	parser.add_argument('-w', '--overwrite', action='store_true', default=False, help='Overwrite directory with default name if it exists')
	# TODO: might need to generalize this to take in multiple input structures?
	parser.add_argument('--start_struct', type=str, default="", help='Starting complex structure to use')
	parser.add_argument('--wt_secstruct', type=str, default="", help='Wildtype secondary structure, if not specified, it will be predicted using RNAfold')
	parser.add_argument('--seq_file', type=str, default="", help='Text file containing sequences to be modeled')
	parser.add_argument('--sfxn', type=str, default="P_overlap_reR-hbond_sc_MOD", help='Rosetta score function to use')
	parser.add_argument('--rosetta_prefix', type=str, default="", help='Prefix for the directory containing rosetta executables')
	parser.add_argument('--relax_cutoff_dist', type=str, default="20.0", help='Relax cutoff distance (for sidechain repacking in mutate_and_score_RNP)')
	parser.add_argument('--protein_pack_reps', type=int, default=10, help='Protein pack reps - for protein mutants (for unbound in mutate_and_score_RNP)')
	parser.add_argument('--Nreps', type=str, default="10", help='Structure pred reps')
	parser.add_argument('--no_min_jumps', action='store_true', default=False, help='Minimize jumps in mutate_and_score_RNP?')
	parser.add_argument('--move_backbone', action='store_true', default=False, help='Move backbone in mutate_and_score_RNP?')
	parser.add_argument('--move_protein_backbone', action='store_true', default=False, help='Move protein backbone in mutate_and_score_RNP?')
	parser.add_argument('--submit_qsub_files', action='store_true', default=False, help='Make and submit qsub files')
	args = parser.parse_args()
	main( args )
