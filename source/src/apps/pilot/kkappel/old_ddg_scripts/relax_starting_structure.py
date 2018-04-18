import os
from sys import exit
import argparse
import subprocess

#Relaxation commands
relax_commands = [
'minimize_with_cst -ignore_zero_occupancy false  -in:file:s ../../START_STRUCT_PREF.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min1  -in:file:fullatom  -ddg::constraint_weight 1.0  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min1.START_STRUCT_PREF_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min2  -in:file:fullatom  -ddg::constraint_weight 0.8  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min2.min1_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min3  -in:file:fullatom  -ddg::constraint_weight 0.6  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min3.min2_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min4  -in:file:fullatom  -ddg::constraint_weight 0.4  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min4.min3_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min5  -in:file:fullatom  -ddg::constraint_weight 0.2  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min5.min4_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min6  -in:file:fullatom  -ddg::constraint_weight 0.1  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min6.min5_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min7  -in:file:fullatom  -ddg::constraint_weight 0.05  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:s min7.min6_0001.pdb  -score:weights SFXN  -ddg:out_pdb_prefix min8  -in:file:fullatom  -ddg::constraint_weight 0.0  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false',

'mutate_and_score_RNP -ignore_zero_occupancy false -s min8.min7_0001.pdb -score:weights SFXN -move_backbone true -move_protein_backbone true -out_prefix min_START_STRUCT_PREF_ -relax_cutoff_dist 200.0 -min_jumps true -protein_pack_reps 0 > relax_1.log',

'mutate_and_score_RNP -ignore_zero_occupancy false -s min_START_STRUCT_PREF_wildtype_bound.pdb -score:weights SFXN -move_backbone true -move_protein_backbone true -out_prefix min_again_START_STRUCT_PREF_ -relax_cutoff_dist 200.0 -min_jumps true -protein_pack_reps 0 > relax_2.log' 
]

def main( args ):
	# First make a directory for doing the relaxation
	run_dir = os.getcwd()
	start_struct_pref = args.start_struct.replace('.pdb','')
	relax_dir = run_dir + '/relax_%s' %(start_struct_pref)
	if not os.path.exists(relax_dir):
		os.mkdir( relax_dir )
	else: 
		print "Directory relax_%s already exists!" %(start_struct_pref)
		exit( 0 )
	os.chdir( relax_dir )
	if not args.submit_qsub_files:
		# write a file that contains all commands, can be run locally to do relaxation
		command_file = open('relax_commands.txt', 'w')
	job_ids = []
	for d in range(1,21):
		os.mkdir( relax_dir + '/' + str(d) )
		if args.submit_qsub_files:
			os.chdir( relax_dir + '/' + str(d) )
			# Make the qsub file
			with open('qsub_initial_min.sh', 'w') as qfil:
				# write the qsub header
				qfil.write('#!/bin/bash\n')
				qfil.write('#PBS -N relax_%d\n' %(d))
				qfil.write('#PBS -o %s/%d/relax_job.out\n' %(relax_dir,d))
				qfil.write('#PBS -e %s/%d/relax_job.err\n' %(relax_dir,d))
				qfil.write('#PBS -m n\n')
				qfil.write('#PBS -M nobody@stanford.edu\n')
				qfil.write('#PBS -l walltime=24:00:00\n\n')
				qfil.write('cd %s/%d\n\n' %(relax_dir,d))
				# write the qsub commands
				for command in relax_commands:
					command_di = command.replace('START_STRUCT_PREF',start_struct_pref)
					command_d = command_di.replace('SFXN',args.sfxn)
					qfil.write( args.rosetta_prefix + command_d + '\n' )
				# Submit the qsub files
			try:
				qsub_output = subprocess.check_output('qsub qsub_initial_min.sh', shell=True)
			except:
				qsub_output = subprocess.check_output('qsub qsub_initial_min.sh', shell=True)
			job_ids.append( qsub_output.split('.')[0] )
		else:
			command_file.write('cd %s/%d\n' %(relax_dir, d))
			for command in relax_commands:
				command_di = command.replace('START_STRUCT_PREF',start_struct_pref)
				command_d = command_di.replace('SFXN',args.sfxn)
				command_file.write( args.rosetta_prefix + command_d + '\n' )
		
	if not args.submit_qsub_files:
		# Find the final lowest scoring structure (in the command_file)
		command_file.write('python ~/bin/find_lowest_scoring_relaxed_struct.py --relax_dir %s --start_struct_pref %s' %(relax_dir, start_struct_pref))
		command_file.close()
	else:
		os.chdir(relax_dir)
		# write one final job that depends on all other jobs that will collect the final lowest scoring structure
		with open('qsub_get_results.sh', 'w') as f:
			f.write('#!/bin/bash\n')
			f.write('#PBS -N relax_results\n')
			f.write('#PBS -o %s/relax_results.out\n' %(relax_dir))
			f.write('#PBS -e %s/relax_results.err\n' %(relax_dir))
			f.write('#PBS -m n\n')
			f.write('#PBS -M nobody@stanford.edu\n')
			f.write('#PBS -l walltime=24:00:00\n\n')
			all_job_ids = ':'.join(job_ids)
			f.write('#PBS -W depend=afterany:%s\n' %(all_job_ids))
			f.write('cd %s\n\n' %(relax_dir))
			# Just run a python script that will get the lowest scoring structure
			# TODO: where to put this python script
			f.write('python ~/bin/find_lowest_scoring_relaxed_struct.py --relax_dir %s --start_struct_pref %s' %(relax_dir, start_struct_pref))
		try:
			output = subprocess.check_output('qsub qsub_get_results.sh', shell=True)
		except:
			output = subprocess.check_output('qsub qsub_get_results.sh', shell=True)

	

if __name__ == '__main__':
	print "\n#################################################"
	print "WARNING: Please use updated versions of these scripts in apps/public/rnp_ddg/"
	print "#################################################\n"
	parser=argparse.ArgumentParser(description="Relax an experimental structure in preparation for ddG predictions")
	parser.add_argument('-s', '--start_struct', type=str, default="", help='Starting structure file')
	parser.add_argument('--submit_qsub_files', action='store_true', default=False, help='Make and submit qsub files')
	parser.add_argument('--rosetta_prefix', type=str, default="", help='Prefix for the directory containing rosetta executables')
	parser.add_argument('--sfxn', type=str, default="P_overlap_reR-hbond_sc_MOD", help='Rosetta score function to use')
	args = parser.parse_args()
	main( args )

