import os
import argparse
import subprocess
import mutant_modeler
from time import sleep

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

	cmd_opts += '-restore_talaris_behavior '
	
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

	local_file = open('ALL_COMMANDS', 'w')
	job_IDs = {}
	# Get command lines for all the mutants
	for mutant_index in mutant_indices:
		local_file_m = open('COMMAND_%d' %(mutant_index), 'w')
		setup_commands, model_cmd_files, final_ddg_commands, depend_indices = modeler.get_command_lines( mutant_index )
		# Make qsub files, if specified	
		run_dir = script_run_dir + '/' + output_dir + '/' + str(mutant_index)

		# Write to a single file that can be run locally
		# Change to the appropriate directory
		local_file.write('cd ' + script_run_dir + '/' + output_dir + '/' + str(mutant_index) + '\n')
		local_file_m.write('cd ' + script_run_dir + '/' + output_dir + '/' + str(mutant_index) + '\n')
		local_file.write(setup_commands + '\n')
		local_file_m.write(setup_commands + '\n')
		local_file.write('cd ' + script_run_dir + '/' + output_dir + '/' + str(mutant_index) + '\n')
		local_file_m.write('cd ' + script_run_dir + '/' + output_dir + '/' + str(mutant_index) + '\n')
		local_file.write(final_ddg_commands + '\n')
		local_file_m.write(final_ddg_commands + '\n')
		
		local_file.write( 'cd ' + script_run_dir + '\n')
		local_file_m.write( 'cd ' + script_run_dir + '\n')
		# Write nearest neighbor energy


		# Figure out the free energy of the sequence (want to write this to outfile)
		ftemp = open('temp_seq.txt', 'w')
		ftemp.write( modeler.mutant_seqs[mutant_index][0] + '\n')
		ftemp.close()
		# Use RNAfold to get the free energy
		p_min = ''
		E_min = ''
		fe = ''
		if args.skip_RNAfold:
			fe = '0.0'
			E_min = '0.0'
			p_min = '0.0'
		else:
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
		local_file_m.close()


	local_file.close()

	print "################################################"
	print "Finished setting up ddG calculations!"
	print "All files are in %s" %(output_dir)
	print "To run type: "
	print "source %s/ALL_COMMANDS" %(output_dir)
	print ""
	print "Alternatively, commands are separated out in the COMMAND_* files."
	print "(This is especially useful if you are running on a cluster: each command can be run separately.)"
	print "These can also be run by typing e.g.: source %s/COMMAND_0" %(output_dir)

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
	parser=argparse.ArgumentParser(description="Set up ddg calculations for list of RNP mutants")
	parser.add_argument('--high_res', action='store_true', default=False, help='Make loop mutations using stepwise, all others threaded')
	parser.add_argument('--med_res', action='store_true', default=False, help='Make non-canonical base pair mutations using farfar, all others threaded')
	parser.add_argument('--low_res', action='store_true', default=False, help='Make all mutations with simple threading protocol (DEFAULT METHOD)')
	parser.add_argument('--skip_RNAfold', action='store_true', default=False, help='Do not do RNAfold calculation.')
	parser.add_argument('-t', '--tag', type=str, default="", help='Tag for naming files and directories')
	parser.add_argument('-w', '--overwrite', action='store_true', default=False, help='Overwrite directory with default name if it exists')
	# TODO: might need to generalize this to take in multiple input structures?
	parser.add_argument('--start_struct', type=str, default="", help='Starting complex structure to use')
	parser.add_argument('--wt_secstruct', type=str, default="", help='Wildtype secondary structure, if not specified, it will be predicted using RNAfold')
	parser.add_argument('--seq_file', type=str, default="", help='Text file containing sequences to be modeled')
	parser.add_argument('--sfxn', type=str, default="rnp_ddg", help='Rosetta score function to use')
	parser.add_argument('--rosetta_prefix', type=str, default="", help='Prefix for the directory containing rosetta executables')
	parser.add_argument('--relax_cutoff_dist', type=str, default="20.0", help='Relax cutoff distance (for sidechain repacking in rnp_ddg)')
	parser.add_argument('--protein_pack_reps', type=int, default=10, help='Protein pack reps - for protein mutants (for unbound in rnp_ddg)')
	parser.add_argument('--Nreps', type=str, default="10", help='Structure pred reps')
	parser.add_argument('--no_min_jumps', action='store_true', default=False, help='Minimize jumps in rnp_ddg?')
	parser.add_argument('--move_backbone', action='store_true', default=False, help='Move backbone in rnp_ddg?')
	parser.add_argument('--move_protein_backbone', action='store_true', default=False, help='Move protein backbone in rnp_ddg?')
	args = parser.parse_args()
	main( args )
