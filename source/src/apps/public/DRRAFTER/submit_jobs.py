#!/usr/bin/env python
import os
import argparse

def setup_out_dir( out_dir, njobs ):
	if not os.path.exists( out_dir ):
		os.mkdir( out_dir )
	for i in range( njobs ):
		if not os.path.exists( out_dir + '/' + str(i) ):
			os.mkdir( out_dir + '/' + str(i) )

def make_job_files_and_submit( command_file, out_dir, out_pref, njobs, template_submission_script, queue_command, no_submit ):
	# get the template submission script lines
	if not os.path.exists( template_submission_script ):
		print( "{file} does not exist!".format(file=template_submission_script) )
		exit()
	template_submission_script_lines = ''
	for line in open( template_submission_script ):
		template_submission_script_lines += line
	
	# make a directory to keep all of the job files
	if not os.path.exists( 'job_files' ):
		os.mkdir( 'job_files' )
	command_file_lines = ''
	for line in open( command_file ):
		if line.startswith( '#' ): continue
		if len( line.replace( '\n', '' ) ) == 0: continue
		base_command = line.split('@')[0]
		command_file_lines += base_command
		if len( line.split( '@' ) ) > 1:
			flags_file = line.split('@')[1].replace('\n','')
			for fline in open( flags_file ):
				command_file_lines += fline.replace( '\n', '') + ' '

	working_dir = os.getcwd()
	all_job_files = []
	for i in range( njobs ):
		submission_script_lines = template_submission_script_lines.replace( 'JOB_NAME', out_dir + '_' + str(i) )
		submission_script_lines += '\ncd ' + working_dir + '\n'
		job_file_name = 'job_files/job{i}.sh'.format( i=i )
		all_job_files.append( job_file_name )
		with open( job_file_name, 'w' ) as f:
			f.write( submission_script_lines )
			out_dir_i = out_dir + '/' + str( i ) + '/'
			command_line_i = command_file_lines.replace( '-out:file:silent ', '-out:file:silent ' + out_dir_i )
			command_line_i = command_line_i.replace( '-silent ', '-out:file:silent ' + out_dir_i )
			f.write( command_line_i )

	# and submit the jobs
	for job_file in all_job_files:
		if no_submit:
			print( '{submit_command} {file}'.format( submit_command=queue_command, file=job_file ) )
		else:
			os.system( '{submit_command} {file}'.format( submit_command=queue_command, file=job_file ) )

def submit_jobs( args ):
	# is this a final round?
	final_round = False
	if "FINAL" in args.curr_round:
		final_round = True
	
	# if it's not a final round, then we need to submit jobs for all possible alignments
	# these fits are in a file: {out_pref}_auto_fits.txt 
	if not os.path.exists( '{out_pref}_auto_fits.txt'.format( out_pref= args.out_pref ) ):
		print( "Cannot find {out_pref}_auto_fits.txt".format( out_pref= args.out_pref ) )
		print( "Are you running from the correct directory?" )
		print( "Or did you set up manually and forget to create this file?" )
		exit()
	all_fits = []
	for line in open( '{out_pref}_auto_fits.txt'.format( out_pref= args.out_pref ) ):
		all_fits.append( line.replace(' ', '' ).replace('\n','') )

	# set up the jobs for each fit
	if final_round:
		command_file = 'command_' + args.out_pref + '_' + args.curr_round
		command_file_half1 = command_file + '_half1'
		command_file_half2 = command_file + '_half2'
		# and might have half maps if it's the final round
		# check if the _half1 and _half2 files exist and use them if they do
		if os.path.exists( command_file_half1 ) and os.path.exists( command_file_half2 ):
			# set up for the halves
			# make the output directory
			out_dir1 = 'out_' + args.out_pref + '_' + args.curr_round + '_half1'
			out_dir2 = 'out_' + args.out_pref + '_' + args.curr_round + '_half2'
			setup_out_dir( out_dir1, args.njobs ) 
			setup_out_dir( out_dir2, args.njobs ) 
			make_job_files_and_submit( command_file_half1, out_dir1, args.out_pref, args.njobs, 
				args.template_submission_script, args.queue_command, args.no_submit )
		else:
			if not os.path.exists( command_file ):
				print( "Cannot find command file %s" %(command_file) )
				exit()
			# set up for the single command file
			out_dir = 'out_' + args.out_pref + '_' + args.curr_round
			setup_out_dir( out_dir, args.njobs )
			make_job_files_and_submit( command_file, out_dir, args.out_pref, args.njobs, 
				args.template_submission_script, args.queue_command, args.no_submit )
	elif "SINGLE_FIT" in all_fits:
		command_file = 'command_' + args.out_pref + '_' + args.curr_round
		out_dir = 'out_' + args.out_pref + '_' + args.curr_round
		setup_out_dir( out_dir, args.njobs )
		make_job_files_and_submit( command_file, out_dir, args.out_pref, args.njobs, 
			args.template_submission_script, args.queue_command, args.no_submit )
	else:
		for fit in all_fits:
			command_file = 'command_' + args.out_pref + '_' + fit + '_' + args.curr_round
			out_dir = 'out_' + args.out_pref + '_' + fit + '_' + args.curr_round
			setup_out_dir( out_dir, args.njobs )
			make_job_files_and_submit( command_file, out_dir, args.out_pref, args.njobs, 
				args.template_submission_script, args.queue_command, args.no_submit )
		

if __name__ == '__main__':
	parser = argparse.ArgumentParser( description='Submit jobs for auto-DRRAFTER' )
	parser.add_argument( '-out_pref', type=str)
	parser.add_argument( '-curr_round', type=str )
	parser.add_argument( '-njobs', type=int )
	parser.add_argument( '-template_submission_script', type=str )
	parser.add_argument( '-queue_command', type=str )
	parser.add_argument( '-no_submit', action='store_true', default=False )
	args = parser.parse_args()
	submit_jobs( args )
