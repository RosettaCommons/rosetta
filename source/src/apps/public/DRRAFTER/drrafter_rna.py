import DRRAFTER_util
import os
import mrcfile
import numpy as np
import glob
import argparse

def get_density_per_residue( model, density_map ):

	density_per_residue = {}

	with mrcfile.open( density_map ) as mrc:
		origin_v1 = mrc.header.origin.x
		origin_v2 = mrc.header.origin.y
		origin_v3 = mrc.header.origin.z
		v1_size = mrc.voxel_size.x
		v2_size = mrc.voxel_size.y
		v3_size = mrc.voxel_size.z
	
		def get_density_at_point( point ):
	
			ind1 = int( point[2]/v1_size + origin_v1)
			ind2 = int( point[1]/v2_size + origin_v2)
			ind3 = int( point[0]/v3_size + origin_v3)
			if ind1 >= mrc.header.nx or ind2 >= mrc.header.ny or ind3 >= mrc.header.nz:
				#print( "point out of map!" )
				density = 0.0
			else:
				density = mrc.data[ ind1, ind2, ind3 ]
				# zero out negative density for now
				if density < 0.0:
					density = 0.0
	
			return density
	
		# go through a bunch of models
		dens_array = mrc.data.copy()
		indices_below_zero = dens_array < 0.0
		dens_array[ indices_below_zero ] = 0.0
		DENS_THR = np.mean( dens_array ) + 1*np.std( dens_array )
		# go through all the atoms in the PDB
		nres = 0
		dens_values = []
		per_res_dens_values = []
		low_dens_residues = []
		last_resnum = 'NA'
		for line in open( model ):
			if not line.startswith( 'ATOM' ): continue
			resnum = line.split()[5]
			if resnum != last_resnum:
				nres += 1
				if last_resnum == 'NA':
					last_resnum = resnum
					residue_dens_values = []
				else:
					mean_residue_dens_val = np.mean( residue_dens_values )
					if mean_residue_dens_val < DENS_THR: low_dens_residues.append( last_resnum )
					per_res_dens_values.append( [last_resnum, mean_residue_dens_val] ) 
					density_per_residue[ last_resnum ] = mean_residue_dens_val
					residue_dens_values = []
					last_resnum = resnum
				
			x = float(line[29:38].replace(' ',''))
			y = float(line[38:46].replace(' ',''))
			z = float(line[46:54].replace(' ',''))
		
			dens = get_density_at_point( [x,y,z] )
			residue_dens_values.append( dens )
			
		# get values for last residue
		mean_residue_dens_val = np.mean( residue_dens_values )
		if mean_residue_dens_val < 0.05: low_dens_residues.append( resnum )
		per_res_dens_values.append( [last_resnum, mean_residue_dens_val] ) 
		density_per_residue[ resnum ] = mean_residue_dens_val
		nres += 1

	return density_per_residue

def write_flags_extend_last_round( outfile_basename, output_tag, last_flags, last_round, nstruct, test ):

	new_flags_file = 'flags_' + outfile_basename + '_' + output_tag
	with open( new_flags_file, 'w') as f:
		for line in open( last_flags ):
			if line.startswith( '-s '):
				f.write( line )
			elif line.startswith( '-initial_structures ' ):
				f.write( line )
			elif line.startswith( '-dock_chunks '):
				f.write( line )
			elif line.startswith( '-out:file:silent '):
				outfile_basename = outfile_basename
				#if last_round:
				#	outfile_basename += '_Rfinal'
				f.write( '-out:file:silent ' + outfile_basename + '_' + output_tag+ '.out\n' )
			elif line.startswith( '-minimize_rounds ') and last_round:
				f.write( '-minimize_rounds 2\n' )
			elif line.startswith( '-nstruct ' ) and nstruct != 0:
				f.write( '-nstruct ' + str(nstruct) + '\n' )
			elif line.startswith( '-dock_into_density' ) and last_round:
				f.write( '-dock_into_density true\n' )
			else:
				f.write( line )
	if test:
		f.write( '-testing:INTEGRATION_TEST\n' )

	return new_flags_file

def write_command_file( command_file_name, flags_file, rosetta_dir, rosetta_ext ):

	with open( command_file_name, 'w') as f:
		f.write( rosetta_dir + 'rna_denovo' + rosetta_ext + '@'+ flags_file )
		f.write( '\n' )


def setup_next_round( models, last_flags, output_tag, outfile_basename, do_not_recalc_convergence=False, 
			last_round=False, skip_clash_check=False, dens_thr=-1.0, nstruct=0, chunk_res=[],
			final_round_2=False, map_to_use='', rosetta_dir='', rosetta_ext='', test=False):

	# figure out the convergence per region
	# define regions based on secondary structure graph
	
	# find the secondary structure from the flags file from the previous round
	for line in open( last_flags ):
		if line.startswith( '-secstruct_file' ):
			secstruct_file = line.split()[1].replace('\n','')
		elif line.startswith( '-edensity:mapfile' ):
			if map_to_use != '':
				dens_map = map_to_use
			else:
				dens_map = os.path.abspath( os.path.dirname( last_flags ) ) + '/' + line.split()[1].replace('\n','')
		elif line.startswith( '-s ' ):
			input_chunks = line.split()[1:]

	for line in open( os.path.abspath( os.path.dirname( last_flags ) ) + '/' + secstruct_file ):
		ss = line.replace('\n','')
		break
	
	G, stem_conn = DRRAFTER_util.graph_from_ss( ss )

	convergence_per_node = {}

	name = outfile_basename

	all_nodes = {}
	# need to modify nodes slightly based on "fixed chunks"
	fixed_num = 0
	all_fixed_resis = []
	if chunk_res:
		for chunk in chunk_res:
			# figure out the residues
			node_name = 'fixed' + str(fixed_num)
			node_resis = [] 
			for resis in chunk.split():
				start_resi = int(resis.split(':')[1].split('-')[0])
				if '-' in resis:
					stop_resi = int(resis.split(':')[1].split('-')[1])
				else:
					stop_resi = int(resis.split(':')[1].split('-')[0])
				for i in range( start_resi, stop_resi+1):
					node_resis.append( i )
					all_fixed_resis.append( i )
			all_nodes[ node_name ] = node_resis
			fixed_num += 1
	for node in G.nodes:
		new_resis_this_node = []
		for r in G.nodes[ node ]['residues']:
			if r not in all_fixed_resis:
				new_resis_this_node.append( r )
		if len( new_resis_this_node ) > 0:
			all_nodes[ node ] = new_resis_this_node
	
	for node in all_nodes.keys():
		# if there's only one or two residues, just skip it:
		if len( all_nodes[ node ] ) <= 2: continue
	
		#print( node, all_nodes[ node ] )
	
		# residues (offset from rosetta #ing by 1)
		residue_str = ''
		for r in all_nodes[ node ]:
			residue_str += str(r+1)
			residue_str += ' '
	
		# figure out the convergence over this region
		model_str = ''
		for x in models:
			model_str += x
			model_str += ' '
		convergence_file = name + '_convergence_node_' + node + '.txt'
		if not (do_not_recalc_convergence and os.path.exists( convergence_file )):
			if test:
				os.system( '%s/drrafter_error_estimation%s -s %s -rmsd_res %s -rmsd_nosuper -rna_rmsd false -mute core basic -testing:INTEGRATION_TEST > %s' %( rosetta_dir, rosetta_ext, model_str, residue_str, convergence_file ) )
			else:
				os.system( '%s/drrafter_error_estimation%s -s %s -rmsd_res %s -rmsd_nosuper -rna_rmsd false -mute core basic > %s' %( rosetta_dir, rosetta_ext, model_str, residue_str, convergence_file ) )
		for line in open( convergence_file ):
			if "convergence" in line:
				conv = float( line.split()[-1] )
				break
		convergence_per_node[ node ] = conv
	
	# go through the best scoring model and figure out which regions are fully sitting in density
	# get the density per residue for a set of models
	dens_per_res_per_model = {}
	for model in models:
		dens_per_res = get_density_per_residue( model, dens_map )
		dens_per_res_per_model[ model ] = dens_per_res
		
	# what's the threshold density value?
	with mrcfile.open( dens_map ) as mrc:
		dens_array = mrc.data.copy()
		indices_below_zero = dens_array < 0.0
		dens_array[ indices_below_zero ] = 0.0
		if dens_thr == -1:
			DENS_THR = 0.4*np.max( dens_array ) 
		else:
			DENS_THR = dens_thr
	
	print("Density threshold: %0.3f" %( DENS_THR))
		
	density_values_per_node_per_model = {}
	density_values_per_node_per_model_min = {}
	for node in all_nodes.keys():
		#if len( G.nodes[ node ]['residues'] ) <= 4: continue
		density_values_per_node_per_model[ node ] = {}
		density_values_per_node_per_model_min[ node ] = {}
	
		# is this whole element in the density map?
		# check for each model
		for model in models:
			density_values_node = []
			for res in all_nodes[ node ]:
				if int( res ) != res:  # to get around a weird feature of converting ss to graph (sometimes 0.5 residues...)
					continue
				dens = dens_per_res_per_model[ model ][ str(res) ]
				density_values_node.append( dens )
		
			if len( density_values_node ) < 1:
				density_values_per_node_per_model[ node ][ model ] = 0.0
			else:
				mean_density_value_node = np.mean( density_values_node )
				min_density_value_node = np.min( density_values_node )
				density_values_per_node_per_model[ node ][ model ] = mean_density_value_node
				density_values_per_node_per_model_min[ node ][ model ] = min_density_value_node

	# go through each node and figure out if it passes the convergence check
	# if yes, then figure out whether the best fit in density passes the density threshold
	CONV_THR = 20.0
	CONV_THR_SMALL_NODES = 7.0
	# different thresholds for different numbers of residues
	CONV_THR_12 = 25.0
	CONV_THR_10 = 20.0
	CONV_THR_8 = 10.0
	CONV_THR_6 = 7.0
	backup_models_best_fit = {}
	# loops need to be really converged in order for us to freeze them out! (will included as a fixed chunk if passes this threshold)
	if last_round:
		CONV_THR_LOOP = 1.5
	else:
		CONV_THR_LOOP = 4.0
	FREEZE_OUT_EVEN_IF_POOR_DENSITY_OR_FEW_RES = 0.8
	best_fit_chunks_and_residues = {}
	for node in all_nodes.keys():
		nres = len( all_nodes[ node ] )
		if len( all_nodes[ node ] ) <= 4:
			if node not in convergence_per_node.keys(): 
				continue
			if convergence_per_node[ node ] > FREEZE_OUT_EVEN_IF_POOR_DENSITY_OR_FEW_RES and not last_round: 
				continue
		# check convergence
		if not (node.startswith( 'S' ) or node.startswith('fixed')):
			if convergence_per_node[ node ] > CONV_THR_LOOP: continue
		elif nres <= 6:
			# smaller helices better be really locked in
			# otherwise we get stuck with things that aren't really correct
			if convergence_per_node[ node ] > CONV_THR_6: continue
		elif nres <= 8:
			if convergence_per_node[ node ] > CONV_THR_8: continue
		elif nres <= 10:
			if convergence_per_node[ node ] > CONV_THR_10: continue
		else:
			if convergence_per_node[ node ] > CONV_THR_12: continue
		if not last_round:
			best_fit = sorted( density_values_per_node_per_model[ node ].items(), key=lambda x: x[1], reverse=True )[0]
		else:
			# take the node from the best scoring model
			for x in density_values_per_node_per_model[ node ].items():
				if "out.1.pdb" in x[0]: # a bit of a hack...
					best_fit = x
					break
		dens_val_best_fit = best_fit[1]
		model_best_fit = best_fit[0]
		# what's the minimum density value per residue (we want every residue to sit in density!)
		min_dens = density_values_per_node_per_model_min[ node ][ model_best_fit ]
		#print node, dens_val_best_fit, min_dens
		if (dens_val_best_fit > DENS_THR and last_round) or (min_dens > DENS_THR and dens_val_best_fit > DENS_THR) or (convergence_per_node[ node ] < FREEZE_OUT_EVEN_IF_POOR_DENSITY_OR_FEW_RES and not last_round):
		#if dens_val_best_fit > DENS_THR or (convergence_per_node[ node ] < FREEZE_OUT_EVEN_IF_POOR_DENSITY_OR_FEW_RES and not last_round):
			#print node, model_best_fit, dens_val_best_fit, "conv", convergence_per_node[ node ], "nres", len( all_nodes[ node ] ), "min dens", min_dens
			# extract this node from the model
			residue_str = ''
			for res in all_nodes[ node ]:
				residue_str += 'A:'
				residue_str +=  str(res)
				residue_str +=  ' '
			#os.system( 'pdbslice.py %s -subset %s %s_region_%s_' %(model_best_fit, residue_str, output_tag, node ) )
			DRRAFTER_util.pdbslice( model_best_fit, "subset", residue_str, '%s_region_%s_' %(output_tag, node))
			best_fit_chunks_and_residues[ output_tag + '_region_' + node + '_' + model_best_fit ] = all_nodes[ node ]
		elif convergence_per_node[ node ] < FREEZE_OUT_EVEN_IF_POOR_DENSITY_OR_FEW_RES:
			#print "BACKUP NODE", node
			residue_str = ''
			for res in all_nodes[ node ]:
				residue_str += 'A:'
				residue_str +=  str(res)
				residue_str +=  ' '
			#os.system( 'pdbslice.py %s -subset %s %s_region_%s_' %(model_best_fit, residue_str, output_tag, node ) )
			DRRAFTER_util.pdbslice( model_best_fit, "subset", residue_str, '%s_region_%s_' %( output_tag, node))
			backup_models_best_fit[ output_tag + '_region_' + node + '_' + model_best_fit ] = all_nodes[ node ]
		#else:
			#print "skipped", node

	if len( best_fit_chunks_and_residues.keys() ) < 1:
		best_fit_chunks_and_residues = backup_models_best_fit

	# figure out if any of the chunks clash with each other -- if they do, omit both (unless one was fit in the previous round -- then keep that one)
	with open( 'tmp_fa_rep.wts' ,'w' ) as f:
		f.write( 'fa_rep 1.0\n' )

	# get residues in the first input chunk, so we know which chunks we need to keep despit potential clashes

	chunk0 = input_chunks[0]
	# figure out what residues are in this chunk
	resnums_in_chunk0 = []
	for line in open( os.path.abspath( os.path.dirname( last_flags ) ) + '/' + chunk0 ):
		if not line.startswith( 'ATOM' ): continue
		resnum = int(line.split()[5])
		resnums_in_chunk0.append( resnum )

	resnums_in_chunk0_set = set( resnums_in_chunk0 )

	CLASH_SCORE_CUTOFF = 100.
	clash_chunks = []
	best_fit_chunks_to_remove = []
	if not skip_clash_check:
		best_fit_chunks_list = list(best_fit_chunks_and_residues.keys())
		for i1 in range( len(best_fit_chunks_list) ):
			for i2 in range( i1+1, len(best_fit_chunks_list) ):
				score = 0.0
				pdb1 = best_fit_chunks_list[ i1 ]
				pdb2 = best_fit_chunks_list[ i2 ]
				os.system( 'cat %s %s > tmp_cat_chunk.pdb' %( pdb1, pdb2 ) )
				#os.system( '~/bin/better_reorder_pdb.py tmp_cat_chunk.pdb > quiet' )
				DRRAFTER_util.better_reorder_pdb( 'tmp_cat_chunk.pdb' )
				# score it with rna_score
				if test:
					os.system( '%s/rna_score%s -s tmp_cat_chunk.REORDER.pdb -out:file:silent tmp.out -overwrite true -score:weights tmp_fa_rep.wts -mute core basic -testing:INTEGRATION_TEST > quiet' %(rosetta_dir, rosetta_ext) )
				else:
					os.system( '%s/rna_score%s -s tmp_cat_chunk.REORDER.pdb -out:file:silent tmp.out -overwrite true -score:weights tmp_fa_rep.wts -mute core basic > quiet' %(rosetta_dir, rosetta_ext) )
				os.system( 'grep "SCORE" tmp.out > tmp.sc' )
				for line in open( 'tmp.sc' ):
					if "description" in line: continue
					score = float(line.split()[1])
				#if score > 1000.: # changing this cutoff on 02.22.19 -- I double checked that it only affects 3DIL R3 and R4
				if score > CLASH_SCORE_CUTOFF:
					# call this a clash! 
					clash_chunks.append( pdb1 )
					clash_chunks.append( pdb2 )
					#print "Chunk", pdb1, " and Chunk", pdb2, "CLASH ( score:", score, ")"
				os.remove( 'tmp_cat_chunk.pdb' )
				os.remove( 'tmp_cat_chunk.REORDER.pdb' )
				os.remove( 'tmp.out' )
				os.remove( 'tmp.sc' )
				os.remove( 'quiet' )
	
		os.remove( 'tmp_fa_rep.wts' )
	

		#### this only works if it's not the last round (then we're not really sure what the first chunk is)
		if not last_round:
			for chunk in clash_chunks:
				chunk_residues = best_fit_chunks_and_residues[ chunk ]
				must_keep_this_chunk = False
				for res in chunk_residues:
					if res in resnums_in_chunk0_set:
						must_keep_this_chunk = True
						break
				if not must_keep_this_chunk:
					best_fit_chunks_to_remove.append( chunk )

	best_fit_chunks_and_residues_final = {}
	for chunk in best_fit_chunks_and_residues.keys():
		if chunk not in best_fit_chunks_to_remove:
			best_fit_chunks_and_residues_final[ chunk ] = best_fit_chunks_and_residues[ chunk ]

	
	# set up a run where each of these helices (or loops) starts in this position, but allowed to dock throughout the run
	# get the rest of the helices -- they are also included as input chunks

	all_chunks_new_run = []
	chunks_not_fit_new_run = []

	# figure out which of the original input chunks we also need
	for i, chunk in enumerate(input_chunks):
		# figure out what residues are in this chunk
		resnums_in_chunk = []
		for line in open( os.path.abspath( os.path.dirname( last_flags ) ) + '/' + chunk ):
			if not line.startswith( 'ATOM' ): continue
			resnum = int(line.split()[5])
			resnums_in_chunk.append( resnum )

		resnums_in_chunk_set = set( resnums_in_chunk )

		this_chunk = ''
		for fit_chunk in best_fit_chunks_and_residues_final.keys():
			# figure out if any of the residues in the fit chunk are in this original input chunk,
			# if so, replace the original input chunk
			for r in best_fit_chunks_and_residues_final[ fit_chunk ]:
				if r in resnums_in_chunk_set:
					this_chunk = fit_chunk
					break
			if this_chunk != '': break
		if this_chunk == '': 
			chunks_not_fit_new_run.append( chunk )
			this_chunk = chunk
		all_chunks_new_run.append( this_chunk )


	# Make -s
	s_string = ''
	for x in all_chunks_new_run:
		s_string += x
		s_string += ' '
	# Make -initial_structures
	# concatenate all best_fit_chunks
	all_best_fit_chunks_str = ''
	for x in best_fit_chunks_and_residues_final.keys():
		all_best_fit_chunks_str += x 
		all_best_fit_chunks_str += ' '
	# concatenate into one PDB
	if len( all_best_fit_chunks_str ) < 1:
		if len( backup_models_best_fit.keys() ) < 1:
			print( "Problem with setting up round. Extending last round." )
			# extend the previous round
			new_flags = write_flags_extend_last_round( outfile_basename=outfile_basename,
				output_tag=output_tag, last_flags=last_flags, last_round=last_round, nstruct=nstruct, test=test )
			# and write a command file
			write_command_file( 'command_' + outfile_basename + '_' + output_tag, new_flags, rosetta_dir, rosetta_ext)
			exit( 1 )
		# add a backup chunk
		all_best_fit_chunks_str += backup_models_best_fit.keys()[0] + ' '
	if len( all_best_fit_chunks_str ) > 0:
		#print( 'cat %s > all_fit_%s_%s.pdb' %(all_best_fit_chunks_str, name, output_tag ) )
		os.system( 'cat %s > all_fit_%s_%s.pdb' %(all_best_fit_chunks_str, name, output_tag ) )
		#os.system( '~/bin/better_reorder_pdb.py all_fit_%s_%s.pdb' %(name, output_tag) )
		DRRAFTER_util.better_reorder_pdb( 'all_fit_%s_%s.pdb' %(name, output_tag) )
	initial_structures_str = 'all_fit_%s_%s.REORDER.pdb ' %( name, output_tag)
	for x in chunks_not_fit_new_run:
		initial_structures_str += x
		initial_structures_str += ' '
	# Make -dock_chunks
	dock_chunks_str = all_best_fit_chunks_str

	if last_round:
		# everything in dock_chunks better also be in s_string
		for x in best_fit_chunks_and_residues_final.keys():
			if x not in all_chunks_new_run:
				s_string += x
				s_string += ' '

	#print( "-s",  s_string )
	#print( "-initial_structures",  initial_structures_str )
	#print( "-dock_chunks",  dock_chunks_str )

	new_flags_file = 'flags_' + outfile_basename + '_' + output_tag
	all_residue_string_for_extra_min_res = ''
	with open( new_flags_file, 'w') as f:
		wrote_s = False
		wrote_initial_structures = False
		wrote_dock_chunks = False
		for line in open( last_flags ):
			if line.startswith( '-s '):
				if not last_round or final_round_2:
					f.write( '-s ' + initial_structures_str + '\n' )
				else:
					f.write( '-s ' + s_string + '\n' )
				wrote_s = True
			elif line.startswith( '-initial_structures ' ):
				if last_round and not final_round_2:
					f.write( '-initial_structures ' + initial_structures_str + '\n' )
				wrote_initial_structures = True
			elif line.startswith( '-dock_chunks '):
				if last_round and not final_round_2:
					f.write( '-dock_chunks ' + dock_chunks_str + '\n' )
				wrote_dock_chunks = True
			elif line.startswith( '-out:file:silent '):
				outfile_basename = outfile_basename
				#if last_round:
				#	outfile_basename += '_Rfinal'
				f.write( '-out:file:silent ' + outfile_basename + '_' + output_tag+ '.out\n' )
			elif line.startswith( '-minimize_rounds ') and last_round:
				f.write( '-minimize_rounds 2\n' )
			elif line.startswith( '-nstruct ' ) and nstruct != 0:
				f.write( '-nstruct ' + str(nstruct) + '\n' )
			elif line.startswith( '-dock_into_density' ) and last_round:
				f.write( '-dock_into_density true\n' )
			elif line.startswith( '-fasta' ):
				# get all the residues in case we're doing extra_min_res
				for fasta_line in open( line.split()[1].replace('\n','') ):
					if fasta_line.startswith( '>' ):
						for part in fasta_line.split():
							if ":" not in part: continue
							all_residue_string_for_extra_min_res += part.replace('\n','') + ' '
				f.write( line )
			elif line.startswith( '-edensity:mapfile' ):
				if map_to_use != '':
					f.write( '-edensity:mapfile ' + map_to_use + '\n' )
				else: f.write( line )
			else:
				f.write( line )
		if final_round_2:
			# extra minimization for all residues in the structure in the final final round
			f.write( '-extra_min_res ' + all_residue_string_for_extra_min_res + '\n' )

		if not wrote_s:
			if last_round and not final_round_2:
				f.write( '-s ' + s_string + '\n' )
			else:
				f.write( '-s ' + initial_structures_str + '\n' )
		if not wrote_initial_structures:
			if last_round and not final_round_2:
				f.write( '-initial_structures ' + initial_structures_str + '\n' )
		if not wrote_dock_chunks:
			if last_round and not final_round_2:
				f.write( '-dock_chunks ' + dock_chunks_str + '\n' )
		if test:
			f.write( '-testing:INTEGRATION_TEST\n' )

	# write a command file
	with open( 'command_' + outfile_basename + '_' + output_tag, 'w') as f:
		f.write( rosetta_dir + '/rna_denovo' + rosetta_ext + ' @' + new_flags_file)
		f.write( '\n' )

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="")
	parser.add_argument('-models', type=str, nargs='+', required=True, help="Best scoring models from the last round of modeling")
	parser.add_argument('-nstruct', type=int, default=0, help="nstruct for the flags file, if not set, it will default to value from last round")
	parser.add_argument('-last_flags', type=str, required=True, help="flags file from last round of modeling")
	parser.add_argument('-output_tag', type=str, help="output_tag e.g. R3")
	parser.add_argument('-map_to_use', type=str, default='', help="use this density map instead of the map from the last flags file")
	parser.add_argument('-outfile_basename', type=str, help="base name for the output silent file e.g. 1GID_auto_0")
	parser.add_argument('-do_not_recalc_convergence', action='store_true', default=False, help="Do not recalculate convergence, if it was already calculated")
	parser.add_argument('-last_round', action='store_true', default=False, help="Set up the last round of modeling? (will set up special refinement)")
	parser.add_argument('-final_round_2', action='store_true', default=False, help="Set up the final-round-2 of modeling? (will set up special refinement)")
	parser.add_argument('-skip_clash_check', action='store_true', default=False, help="Don't check for clashes")
	parser.add_argument('-chunk_res', type=str, nargs='+', default=False, help="Don't check for clashes")
	parser.add_argument('-dens_thr', type=float, default=-1.0, help="Threshold average density value for chunk to be fixed -- useful to play around with this if the fixed helices don't look like they fit in the map well")
	parser.add_argument( '-rosetta_directory', required=True, help="Path to Rosetta executables" )
	args = parser.parse_args()
	setup_next_round( models=args.models, nstruct=args.nstruct, last_flags=args.last_flags,
			output_tag=args.output_tag, outfile_basename=args.outfile_basename, 
			do_not_recalc_convergence=args.do_not_recalc_convergence, last_round=args.last_round,
			skip_clash_check=args.skip_clash_check, chunk_res=args.chunk_res, dens_thr=args.dens_thr,
			final_round_2=args.final_round_2, map_to_use=args.map_to_use, rosetta_dir=args.rosetta_directory )
