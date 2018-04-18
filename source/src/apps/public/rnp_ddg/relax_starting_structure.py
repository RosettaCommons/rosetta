import os
from sys import exit
import argparse
import subprocess

#Relaxation commands
relax_commands = [
'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_1.txt  -score:weights SFXN  -ddg:out_pdb_prefix min1  -in:file:fullatom  -ddg::constraint_weight 1.0  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_2.txt  -score:weights SFXN  -ddg:out_pdb_prefix min2  -in:file:fullatom  -ddg::constraint_weight 0.8  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_3.txt  -score:weights SFXN  -ddg:out_pdb_prefix min3  -in:file:fullatom  -ddg::constraint_weight 0.6  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_4.txt  -score:weights SFXN  -ddg:out_pdb_prefix min4  -in:file:fullatom  -ddg::constraint_weight 0.4  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_5.txt  -score:weights SFXN  -ddg:out_pdb_prefix min5  -in:file:fullatom  -ddg::constraint_weight 0.2  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_6.txt  -score:weights SFXN  -ddg:out_pdb_prefix min6  -in:file:fullatom  -ddg::constraint_weight 0.1  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_7.txt  -score:weights SFXN  -ddg:out_pdb_prefix min7  -in:file:fullatom  -ddg::constraint_weight 0.05  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'minimize_with_cst -ignore_zero_occupancy false  -in:file:l ../input_lists/input_8.txt  -score:weights SFXN  -ddg:out_pdb_prefix min8  -in:file:fullatom  -ddg::constraint_weight 0.0  -ddg::cst_dist_cutoff 15.0   -ddg::sc_min_only false -restore_talaris_behavior',

'rnp_ddg -ignore_zero_occupancy false -s min8.min7_0001.pdb -score:weights SFXN -move_backbone true -move_protein_backbone true -out_prefix min_START_STRUCT_PREF_ -relax_cutoff_dist 200.0 -min_jumps true -protein_pack_reps 0 -restore_talaris_behavior > relax_1.log',

'rnp_ddg -ignore_zero_occupancy false -s min_START_STRUCT_PREF_wildtype_bound.pdb -score:weights SFXN -move_backbone true -move_protein_backbone true -out_prefix min_again_START_STRUCT_PREF_ -relax_cutoff_dist 200.0 -min_jumps true -protein_pack_reps 0 -restore_talaris_behavior > relax_2.log' 
]

def main( args ):

	if args.nstructs < 10:
		print "Warning: the recommend number of structures (nstructs) is 100!"

	# First make a directory for doing the relaxation
	run_dir = os.getcwd()
	start_struct_pref = args.start_struct.replace('.pdb','')
	relax_dir = run_dir + '/relax_%s' %(start_struct_pref)
	if not os.path.exists(relax_dir):
		os.mkdir( relax_dir )
	else: 
		print "Directory relax_%s already exists!" %(start_struct_pref)
		exit( 0 )
	command_file = open('ALL_RELAX_COMMANDS', 'w')
	command_file.write( 'cd %s\n' %(relax_dir))
	#os.chdir( relax_dir )
	# write a file that contains all commands, can be run locally to do relaxation
	job_ids = []
	for d in range(1,args.nstructs + 1):
		os.mkdir( relax_dir + '/' + str(d) )
		command_file_d = open('RELAX_COMMAND_%d' %(d), 'w')
		command_file.write('cd %s/%d\n' %(relax_dir, d))
		command_file_d.write('cd %s/%d\n' %(relax_dir, d))
		prefix = args.rosetta_prefix
		if len( prefix ) > 0 and not prefix.endswith('/'):
			prefix += '/'
		for command in relax_commands:
			command_di = command.replace('START_STRUCT_PREF',start_struct_pref)
			command_d = command_di.replace('SFXN',args.sfxn)
			command_file.write( prefix + command_d + '\n' )
			command_file_d.write( prefix + command_d + '\n' )
		command_file_d.close()
	# write the list files
	os.mkdir( '%s/input_lists' %(relax_dir) )
	with open( '%s/input_lists/input_1.txt' %(relax_dir), 'w' ) as f:
		f.write( '../../%s.pdb\n' %(start_struct_pref) )

	with open( '%s/input_lists/input_2.txt' %(relax_dir), 'w') as f:
		f.write( 'min1.%s_0001.pdb\n' %(start_struct_pref) )
		
	for i in range( 3, 9):
		with open( '%s/input_lists/input_%d.txt' %(relax_dir,i), 'w') as f:
			f.write( 'min%d.min%d_0001.pdb\n' %(i-1, i-2) )
		
	# Find the final lowest scoring structure (in the command_file)
	command_file.write( 'cd %s\n' %(run_dir))
	command_file.close()

	print "################################################"
	print "Done setting up relaxation runs!"
	print "See instructions below to run."
	print "################################################"
	print "Commands for running relaxation have been written to ALL_RELAX_COMMANDS"
	print "To run type: "
	print "source ALL_RELAX_COMMANDS"
	print ""
	print "Alternatively, commands are separated out in the RELAX_COMMAND_* files."
	print "(This is especially useful if you are running on a cluster: each command can be run separately.)"
	print "These can also be run by typing e.g.: source RELAX_COMMAND_1"
	

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Relax an experimental structure in preparation for ddG predictions")
	parser.add_argument('-s', '--start_struct', type=str, default="", help='Starting structure file')
	parser.add_argument('--rosetta_prefix', type=str, default="", help='Prefix for the directory containing rosetta executables')
	parser.add_argument('--sfxn', type=str, default="rnp_ddg", help='Rosetta score function to use')
	parser.add_argument('--nstructs', type=int, default=100, help='Number of relaxations to perform (100 is recommended)')
	args = parser.parse_args()
	main( args )
