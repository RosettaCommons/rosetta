#!/usr/bin/env python3
import glob
import os
import argparse
import drrafter_rna
import DRRAFTER_util

def figure_out_next_round( args ):

	# figure out the overall best scoring models from the last round
	# determine whether those models were from final-round-1 
	# (or final-round-2 -- that shouldn't happen, b/c final-round-2 is the very last round)

	this_round_num = int( args.curr_round.replace('FINAL_R','').replace('R','') )

	model_final = args.out_pref + '_all_models_all_fits_FINAL_R' + str( this_round_num )
	model_final += '.out.1.pdb'

	last_round_was_final_R1 = False
	if os.path.exists( model_final ):
		last_round_was_final_R1 = True

	models_str = ''
	for i in range( 1, 11 ):
		model = args.out_pref + '_all_models_all_fits_' + args.curr_round
		model += '.out.' + str(i) + '.pdb'
		if not os.path.exists( model ):
			print( "%s does not exist!" %(model) )
			exit( 1 )
		models_str += model + ' '

	convergence_file = 'convergence_' + args.out_pref + '_all_models_' + args.curr_round + '.txt'
	if args.test:
		os.system( '{rosetta_dir}/drrafter_error_estimation{ext} -s {models_str} -rmsd_nosuper -mute core basic -testing:INTEGRATION_TEST > {fil}'.format( rosetta_dir=args.rosetta_directory, ext=args.rosetta_extension, models_str=models_str, fil=convergence_file ) )
	else:
		os.system( '{rosetta_dir}/drrafter_error_estimation{ext} -s {models_str} -rmsd_nosuper -mute core basic > {fil}'.format( rosetta_dir=args.rosetta_directory, ext=args.rosetta_extension, models_str=models_str, fil=convergence_file ) )

	# get the convergence of the overall best scoring models from the last round
	# write this out to a file
	FIRST_FINAL_ROUND = False
	SECOND_FINAL_ROUND = False
	conv = ''
	for line in open( convergence_file ):
		if "convergence" in line:
			conv = float( line.split()[-1] )
			break
	if conv == '':
		print( "Trouble calculating convergence!" )
		exit( 1 )

	print( "Overall convergence %0.3f" %(conv) )

	# figure out whether the convergence is good enough to go to the "final round"
	if last_round_was_final_R1:
		SECOND_FINAL_ROUND = True
	elif conv < args.convergence_threshold:
		FIRST_FINAL_ROUND = True
	else:
		FIRST_FINAL_ROUND = False

	SKIP_TO_FINAL_ROUND = False
	# determine whether the convergence is decreasing an acceptable amount
	# if not, then switch to first final round, BUT w/ each alignment present in top 10 scoring...
	if this_round_num > 4 and not FIRST_FINAL_ROUND and not SECOND_FINAL_ROUND and not args.do_not_check_convergence_decrease:
		# get the convergence from the last few rounds
		# check that the average change in convergence is > some cutoff?
		# but if the convergence is greater than some number, then even if we're stuck it might not be worth going to the final round?
		convergence_over_last_rounds = []
		for r_num in range( this_round_num - 3, this_round_num + 1):
			for line in open( 'convergence_' + args.out_pref + '_all_models_R' + str(r_num) + '.txt' ):
				if "Mean pairwise RMSD" in line:
					conv = float(line.split()[5].replace('\n',''))
					convergence_over_last_rounds.append( conv )
					break
		# figure out the average change in convergence from round to round
		delta_convergence = []
		delta_convergence_sum = 0.
		for i in range( len(convergence_over_last_rounds) -1 ):
			delta = convergence_over_last_rounds[i+1] - convergence_over_last_rounds[i]
			delta_convergence.append( delta )
			delta_convergence_sum += delta
		#print "delta convergence"
		#print delta_convergence
		avg_delta_convergence = delta_convergence_sum / (len( convergence_over_last_rounds ) - 1)
		#print "average delta convergence"
		#print avg_delta_convergence
		#print "convergence_over_last_rounds"
		#print convergence_over_last_rounds
		if avg_delta_convergence > -1.5:
			# then we're stuck!
			if convergence_over_last_rounds[-1] > 30.:
				# then we're nowhere close, so it's not worth doing a final final round?
				print( "Convergence doesn't seem to be getting better, and we're also not ready to go to the final round, stopping. To override, use option -do_not_check_convergence_decrease" )
			SKIP_TO_FINAL_ROUND = True
			print( "Convergence doesn't seem to be getting better, skipping to final round" )

	# figure out the different fits
	fit_file = '%s_auto_fits.txt' %(args.out_pref)
	if not os.path.exists( fit_file ):
		print( "Can't find the file listing fits: %s" %(fit_file) )
	fits = []
	for line in open( fit_file ):
		fits.append( line.split()[0].replace( '\n', '' ) )

	# for each fit, do the setup
	if SECOND_FINAL_ROUND:
		# if map halves are provided, then we want to model separately into each
		# stuff to figure out
		last_flags = 'flags_' + args.out_pref + '_' + args.curr_round
		models = []
		for i in range(1,11):
			model = args.out_pref + '_all_models_all_fits_' + args.curr_round
			model += '.out.' + str(i) + '.pdb'
			models.append( model )
		use_map_halves = False
		if args.map_half_1 != '' and args.map_half_2 != '':
			use_map_halves = True
		if use_map_halves:
			for i, map_half in enumerate([args.map_half_1, args.map_half_2]):
				output_tag = 'FINAL_R' + str( this_round_num + 1 ) + '_half' + str(i+1)
				drrafter_rna.setup_next_round( last_flags=last_flags, models=models, output_tag=output_tag,
					outfile_basename=args.out_pref, last_round=True, chunk_res=args.chunk_res,
					do_not_recalc_convergence=args.do_not_recalc_convergence, final_round_2=True, map_to_use=map_half,
					rosetta_ext=args.rosetta_extension, rosetta_dir=args.rosetta_directory, test=args.test, dens_thr=args.dens_thr )
		else:
			output_tag = 'FINAL_R' + str( this_round_num + 1 )
			drrafter_rna.setup_next_round( last_flags=last_flags, models=models, output_tag=output_tag,
				outfile_basename=args.out_pref, last_round=True, test=args.test, chunk_res=args.chunk_res,
				do_not_recalc_convergence=args.do_not_recalc_convergence, final_round_2=True,
				rosetta_ext=args.rosetta_extension, rosetta_dir=args.rosetta_directory, dens_thr=args.dens_thr )
	elif FIRST_FINAL_ROUND or SKIP_TO_FINAL_ROUND:
		output_tag = 'FINAL_R' + str( this_round_num + 1 )
		# we need to figure out the best fit
		best_fit = ''
		overall_min_score = 0.
		if "SINGLE_FIT" not in fits:
			for fit in fits:
			        score_file = args.out_pref + '_' + fit + '_' + args.curr_round + '_CAT_ALL_ROUNDS.sc'
			        if not os.path.exists( score_file ):
			                print( "Cannot find score file %s" %( score_file ) )
			        scores = []
			        for line in open( score_file ):
			                if "description" in  line: continue
			                scores.append( float(line.split()[1]) )
			        if  best_fit == '' or min(scores) <  overall_min_score:
			                best_fit = fit
			                overall_min_score = min( scores )
		models = []
		for i in range(1,11):
			model = args.out_pref + '_all_models_all_fits_' + args.curr_round
			model += '.out.' + str(i) + '.pdb'
			models.append( model )
		if "SINGLE_FIT" not in fits:
			last_flags = 'flags_' + args.out_pref + '_' + best_fit + '_' + args.curr_round
		else:
			last_flags = 'flags_' + args.out_pref + '_' + args.curr_round
		drrafter_rna.setup_next_round( last_flags=last_flags, models=models, output_tag=output_tag,
			outfile_basename=args.out_pref, last_round=True, test=args.test,
			do_not_recalc_convergence=args.do_not_recalc_convergence, chunk_res=args.chunk_res,
			rosetta_ext=args.rosetta_extension, rosetta_dir=args.rosetta_directory, dens_thr=args.dens_thr )
	else:
		# need to do the setup for each fit
		output_tag = 'R' + str( this_round_num + 1 )
		for fit in fits:
			if fit=="SINGLE_FIT":
				last_flags = 'flags_' + args.out_pref + '_' + args.curr_round
			else:
				last_flags = 'flags_' + args.out_pref + '_' + fit + '_' + args.curr_round
			models = []
			for i in range(1,11):
				if fit=="SINGLE_FIT":
					model = args.out_pref + '_' + args.curr_round + '_CAT_ALL_ROUNDS.out.%d.pdb' %(i)
				else:
					model = args.out_pref + '_' + fit + '_' + args.curr_round + '_CAT_ALL_ROUNDS.out.%d.pdb' %(i)
				models.append( model )
			if fit=="SINGLE_FIT":
				drrafter_rna.setup_next_round( last_flags=last_flags, models=models, output_tag= output_tag,
					outfile_basename=args.out_pref, do_not_recalc_convergence=args.do_not_recalc_convergence,
					rosetta_ext=args.rosetta_extension, rosetta_dir=args.rosetta_directory, test=args.test,
					chunk_res=args.chunk_res, dens_thr=args.dens_thr )
			else:
				drrafter_rna.setup_next_round( last_flags=last_flags, models=models, output_tag= output_tag,
					outfile_basename=args.out_pref + '_' + fit, do_not_recalc_convergence=args.do_not_recalc_convergence,
					rosetta_ext=args.rosetta_extension, rosetta_dir=args.rosetta_directory, test=args.test,
					chunk_res=args.chunk_res, dens_thr=args.dens_thr )
	
def get_top_ten_models_from_silent_files( silent_files, output_name ):
	# get the top ten models from a list of silent files
	scores_and_names_and_files = []
	for sf in silent_files:
		for line in open( sf ):
			if "description" in line: continue
			if not "SCORE" in line: continue
			score = float(line.split()[1])
			name = line.split()[-1]
			scores_and_names_and_files.append( [score, sf, name] )
	# sort
	top_ten = sorted( scores_and_names_and_files, key= lambda x: x[0] )[0:10]
	for index, model in enumerate(top_ten):
		if args.test:
			os.system( '{rosetta_dir}/extract_pdbs{ext} -load_PDB_components -in:file:silent {silent_file} -in:file:silent_struct_type binary -in:file:fullatom -in:file:tags {name} -mute core basic extract_pdbs -testing:INTEGRATION_TEST'.format(silent_file=model[1], name=model[2], rosetta_dir=args.rosetta_directory, ext=args.rosetta_extension ) )
		else:
			os.system( '{rosetta_dir}/extract_pdbs{ext} -load_PDB_components -in:file:silent {silent_file} -in:file:silent_struct_type binary -in:file:fullatom -in:file:tags {name} -mute core basic extract_pdbs'.format(silent_file=model[1], name=model[2], rosetta_dir=args.rosetta_directory, ext=args.rosetta_extension ) )
		os.system( 'mv {name}.pdb {output_name}{num}.pdb'.format( name=model[2], output_name=output_name, num=index+1 ) )
	# make a score file
	full_score_file = output_name.replace('.out.','') + '.sc'
	for sf in silent_files:
		os.system( 'grep "SCORE" %s > %s' %(sf, full_score_file ) )

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_silent_file_subset( silent_file, subset ):

	#read the silent file
	header = ''
	scores_data = []
	new_score_data = None
	is_header_over = False
	for line in open( silent_file ):
	    if not is_header_over:
	        if len(line) > 5 and line[:5] == 'SCORE':
	            if is_number(line.split()[1]):
	                is_header_over = True
	    if not is_header_over:
	        header += line
	    else:
	        if len(line) > 5 and line[:5] == 'SCORE':
	            if new_score_data is not None:
	                scores_data.append(new_score_data)
	            new_score_data = line
	        else:
	            new_score_data += line
	else:
	    scores_data.append(new_score_data)
	
	if len( scores_data ) < subset:
		print( "ERROR: Not enough models!" )
		print( "You currently have %d, but are requesting %d this round" %(len(scores_data), subset) )
		exit()

	# output the subset
	with open( silent_file.replace( '.out', '_subset.out'), 'w') as outf:
		outf.write(header)
		for index in range( subset ):
			outf.write( scores_data[index] )

def get_results( args ):
	# is this a final round?
	final_round = False
	if "FINAL" in args.curr_round:
		final_round = True
	final_round2 = False
	
	if not os.path.exists( '{out_pref}_auto_fits.txt'.format( out_pref= args.out_pref ) ):
		print( "Cannot find {out_pref}_auto_fits.txt".format( out_pref= args.out_pref ) )
		print( "Are you running from the correct directory?" )
		print( "Or did you set up manually and forget to create this file?" )
		exit()
	all_fits = []
	for line in open( '{out_pref}_auto_fits.txt'.format( out_pref= args.out_pref ) ):
		all_fits.append( line.replace(' ', '' ).replace('\n','') )

	if final_round:
		out_dir = 'out_' + args.out_pref + '_' + args.curr_round
		out_dir_half1 = out_dir + '_half1'
		out_dir_half2 = out_dir + '_half2'
		this_round_num = int(args.curr_round.replace( 'FINAL_R', '' ))
		# and might have half maps if it's the final round
		# check if the _half1 and _half2 files exist and use them if they do
		if os.path.exists( out_dir_half1 ) and os.path.exists( out_dir_half2 ):
			final_round2 = True
			# get results for final round built into half maps
			for half in ['1','2']:
				DRRAFTER_util.easy_cat( '{out_dir}_half{half}'.format( out_dir=out_dir, half=half) )
				#os.system( 'easy_cat.py {out_dir}_half{half}'.format( out_dir=out_dir, half=half) )
				silent_file = args.out_pref + '_' + args.curr_round + '_half' + half + '.out' 
				silent_file_subset = args.out_pref + '_' + args.curr_round + '_half' + half + '_subset.out' 
				get_silent_file_subset( silent_file, args.nmodels )
				silent_files_all_rounds = [ silent_file ]
				output_name = args.out_pref + '_all_models_all_fits_' + args.curr_round + '_half' + half + '.out.'
				get_top_ten_models_from_silent_files( silent_files_all_rounds, output_name )
		else:
			if not os.path.exists( out_dir ):
				print( "Cannot find output directory %s" %(out_dir) )
				exit()
			# get results for the final round, not built into half maps
			DRRAFTER_util.easy_cat( out_dir )
			#os.system( 'easy_cat.py {out_dir}'.format( out_dir=out_dir) )
			silent_file = args.out_pref + '_' + args.curr_round + '.out'
			silent_file_subset = args.out_pref + '_' + args.curr_round + '_subset.out'
			get_silent_file_subset( silent_file, args.nmodels )
			silent_files_all_rounds = []
			# not sure if this is final round 1 or final round 2
			for r in range( 1, this_round_num-1 ):
				for fit in all_fits:
					if fit == "SINGLE_FIT":
						silent_file_r = args.out_pref + '_R' + str(r) + '_subset.out'
					else:
						silent_file_r = args.out_pref + '_' + str(fit) + '_R' + str(r) + '_subset.out'
					silent_files_all_rounds.append( silent_file_r )
			silent_files_all_rounds.append( silent_file_subset )
			silent_file_prev_round_final = args.out_pref + '_FINAL_R' + str(this_round_num-1) + '_subset.out'
			output_name = args.out_pref + '_all_models_all_fits_' + args.curr_round + '.out.'
			if os.path.exists( silent_file_prev_round_final ):
				silent_files_all_rounds.append( silent_file_prev_round_final )
				final_round2 = True
			else:
				for fit in all_fits:
					if fit == "SINGLE_FIT":
						silent_file_prev_round = args.out_pref + '_R' + str(this_round_num-1) + '_subset.out'
					else:
						silent_file_prev_round = args.out_pref + '_' + str(fit) + '_R' + str(this_round_num-1) + '_subset.out'
					silent_files_all_rounds.append( silent_file_prev_round )
			get_top_ten_models_from_silent_files( silent_files_all_rounds, output_name )

	elif "SINGLE_FIT" in all_fits:
		out_dir = 'out_' + args.out_pref + '_' + args.curr_round
		# get the models from this directory
		DRRAFTER_util.easy_cat( out_dir )
		#os.system( 'easy_cat.py {out_dir}'.format( out_dir=out_dir) )
		silent_file = args.out_pref + '_' + args.curr_round + '.out'
		silent_file_subset = args.out_pref + '_' + args.curr_round + '_subset.out'
		get_silent_file_subset( silent_file, args.nmodels )
		silent_files_all_rounds = []
		this_round_num = int(args.curr_round.replace( 'R', '' ))
		for r in range( 1, this_round_num):
			silent_file_r = args.out_pref + '_R' + str(r) + '_subset.out'
			silent_files_all_rounds.append( silent_file_r )
		silent_files_all_rounds.append( silent_file_subset )
		output_name = args.out_pref + '_' + args.curr_round + '_CAT_ALL_ROUNDS.out.'
		get_top_ten_models_from_silent_files( silent_files_all_rounds, output_name )
		# these are also the overall top scoring models b/c there's just a single fit
		output_name = args.out_pref + '_all_models_all_fits_' + args.curr_round + '.out.'
		get_top_ten_models_from_silent_files( silent_files_all_rounds, output_name )
	else:
		silent_files_all_rounds_all_fits = []
		for fit in all_fits:
			out_dir = 'out_' + args.out_pref + '_' + fit + '_' + args.curr_round
			# get the models from this directory
			DRRAFTER_util.easy_cat( out_dir )
			#os.system( 'easy_cat.py {out_dir}'.format( out_dir=out_dir ) )
			silent_file = args.out_pref + '_' + fit + '_' + args.curr_round + '.out'
			silent_file_subset = args.out_pref + '_' + fit + '_' + args.curr_round + '_subset.out'
			# get a subset of models (args.nmodels)
			get_silent_file_subset( silent_file, args.nmodels )
			# get the top scoring models overall for this fit
			silent_files_all_rounds = []
			this_round_num = int(args.curr_round.replace( 'R', '' ))
			for r in range( 1, this_round_num):
				silent_file_r = args.out_pref + '_' + fit + '_R' + str(r) + '_subset.out'
				silent_files_all_rounds.append( silent_file_r )
				silent_files_all_rounds_all_fits.append( silent_file_r )
			silent_files_all_rounds.append( silent_file_subset )
			silent_files_all_rounds_all_fits.append( silent_file_subset )
			output_name = args.out_pref + '_' + fit + '_' + args.curr_round + '_CAT_ALL_ROUNDS.out.'
			get_top_ten_models_from_silent_files( silent_files_all_rounds, output_name )
		output_name = args.out_pref + '_all_models_all_fits_' + args.curr_round + '.out.'
		get_top_ten_models_from_silent_files( silent_files_all_rounds_all_fits, output_name )

	setup_next_round = True
	if final_round2: setup_next_round = False

	return setup_next_round

def get_results_and_setup_next_round( args ):
	setup_next_round = get_results( args )
	if setup_next_round:
		print( "Setting up next round" )
		figure_out_next_round( args )
	else:
		print( "DONE building models for " + args.out_pref )
		# finalize models here (fix numbering --
		# really should keep correct numbering throughout, but this should work for now )
		return

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="")
	parser.add_argument( '-curr_round', type=str )
	parser.add_argument( '-nmodels', type=int, default=2000 )
	parser.add_argument('-out_pref', type=str, help="base name for the output files e.g. 1GID_auto")
	parser.add_argument('-chunk_res', type=str, nargs='+', default=False, help="Don't check for clashes")
	parser.add_argument('-convergence_threshold', type=float, default=10.0, help="Convergence threshold for going to the next round")
	parser.add_argument('-map_half_1', type=str, default='', help="Half map 1, if provided with half map 2, then it will be used in final_round_2 modeling")
	parser.add_argument('-map_half_2', type=str, default='', help="Half map 2, if provided with half map 1, then it will be used in final_round_2 modeling")
	parser.add_argument('-do_not_recalc_convergence', action='store_true', default=False, help="Do not recalculate convergence per region, really should only be used for testing")
	parser.add_argument('-test', action='store_true', default=False, help="For running tests on the testing server")
	parser.add_argument('-do_not_check_convergence_decrease', action='store_true', default=False, help="Do not check convergence decrease (will normally stop running if convergence stays the same over many rounds).")
	parser.add_argument( '-rosetta_directory', required=True, help="Path to Rosetta executables" )
	parser.add_argument( '-rosetta_extension', default='', help="Extension for Rosetta executables e.g. '.linuxgccrelease'" )
	parser.add_argument('-dens_thr', type=float, default=-1.0, help="Threshold average density value for chunk to be fixed -- useful to play around with this if the fixed helices don't look like they fit in the map well")
	args = parser.parse_args()
	get_results_and_setup_next_round( args )
