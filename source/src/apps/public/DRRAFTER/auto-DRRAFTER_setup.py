#!/usr/bin/env python
import numpy as np
import networkx as nx
import os
import itertools
import argparse
import datetime
import mrcfile
import DRRAFTER_util

##################################
# Run variables
##################################
END_NODE_ANGLE_CUTOFF = 0.87

FLAGS_TEMPLATE = '''\
-fasta {fasta}
-secstruct_file {secstruct}
-s {cat_helix} {other_helices}
-ft_close_chains false
-edensity:mapfile {map_file}
-edensity:mapreso {map_reso}
-new_fold_tree_initializer true
-bps_moves false
-out:file:silent {out_pref}_R1.out
-minimize_rna true
-minimize_protein_sc true
-rna_protein_docking true
-rnp_min_first true
-rnp_pack_first true
-cycles {cycles}
-rnp_high_res_cycles 2
-minimize_rounds 1
-nstruct {nstruct}
-dock_into_density false
-ignore_zero_occupancy false
-convert_protein_CEN false
-FA_low_res_rnp_scoring true
-ramp_rnp_vdw true
-docking_move_size 0.5 
-dock_each_chunk_per_chain false
-mute protocols.moves.RigidBodyMover
-mute protocols.rna.denovo.movers.RNA_HelixMover
-use_legacy_job_distributor true
-no_filters
-set_weights linear_chainbreak 20.0
-jump_library_file RNA18_HUB_2.154_2.5.jump 
-vall_torsions RNA18_HUB_2.154_2.5.torsions
-score:weights stepwise/rna/rna_res_level_energy4.wts 
-restore_talaris_behavior
-edensity:cryoem_scatterers
'''

def get_num_neighbors( tree, node ):
	num_neighbors = 0
	for k in tree.neighbors( node ):
		num_neighbors += 1

	return num_neighbors

def get_neighbors( tree, node ):
	neighbors = []
	for k in tree.neighbors( node ):
		neighbors.append( k )

	return neighbors
		
def get_distance( p1, p2 ):
	x = p1[0] - p2[0]
	y = p1[1] - p2[1]
	z = p1[2] - p2[2]
	dist = np.sqrt( x**2 + y**2 + z**2 )
	return dist

def get_angle_btwn( v1, v2 ):
	return np.arccos( np.dot( v1, v2)/ ( np.linalg.norm( v1 ) * np.linalg.norm( v2 ) ) )


def get_density_btwn_points( dens_map, point1, point2, npoints ):

	with mrcfile.open( dens_map ) as mrc:
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
			density = mrc.data[ ind1, ind2, ind3 ]

			return density
	
		density_point1 = get_density_at_point( point1 )
		density_point2 = get_density_at_point( point2 )

		# find the amount of density along a line between 2 points
		density_btwn = []
		# points along the line
		for x_i in range(1,npoints+1):
			pt_btwn1 = x_i*(np.array(point2) - np.array(point1))/npoints + np.array(point1)
			density = get_density_at_point( pt_btwn1 )
			density_btwn.append( density )
	
		avg_dens_btwn = np.mean( density_btwn )
	
	return avg_dens_btwn

def auto_setup_helices( args ):

	# Check that all executables exist
	# we need e2proc3d.py, e2segment3d.py, extract_lowscore_decoys.py

	SEQ = ''
	for line in open( args.fasta ):
		if line.startswith( '>' ): continue
		SEQ += line.replace('\n','').replace(' ','' )
	if not SEQ.islower():
		print( "ERROR: RNA sequence must be specified with lowercase letters (a, c, g, and u).")
		exit()
	NRES = len( SEQ )
	SS = ''
	for line in open( args.secstruct ):
		SS += line.replace('\n','' ).replace(' ', '' )
	NINIT_POINTS = NRES/6
	DISTANCE_CONNECT = 20.0 # hardcoded
	DENS_MAP = '{OUT_PREF}_lp20.mrc'.format( OUT_PREF=args.out_pref )

	# check that the length of the sequence and secondary structure are the same
	if len( SEQ ) != len( SS ):
		print( "ERROR: Length of sequence does not match length of secondary structure!" )
		print( "Secondary structure: %s" %(SS) )
		print( "Sequence: %s" %(SEQ) )
		print( "Make sure that your sequence and secondary structure files do not contain any extra information" )
		print( "(e.g. your sequence file should not contain the secondary structure and the secondary structure file should not contain the sequence)" )
		exit()

	# figure out if there's a stretch of single stranded residues that's greater than or equal to 8 nucleotides:
	USE_ONLY_ONE_END_NODE = args.fit_only_one_helix
	#USE_ONLY_ONE_END_NODE = False
	if SS.count( '........' ) > 0:
		USE_ONLY_ONE_END_NODE = True
		print( "Using only one end node!" )
	
	###################################
	# Step 0: low-pass filter the full map
	###################################
	print( "Low-pass filtering the map to 20A." )
	if not args.test:
		os.system('e2proc3d.py {FULL_DENS_MAP} {DENS_MAP} --process=filter.lowpass.tophat:cutoff_freq=0.05'.format(DENS_MAP=DENS_MAP, FULL_DENS_MAP=args.full_dens_map ))
	if args.just_low_pass: 
		return

	# write a settings file 
	with open( 'settings_{out_pref}.txt'.format( out_pref=args.out_pref), 'w' ) as f:
		# write the date and time
		if not args.test:
			f.write( 'Auto helix setup run at ' + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + '\n' )
		f.write( 'MAP_THR: %0.3f\n' %( args.map_thr ) )
		f.write( 'DENS_MAP: %s\n' %(DENS_MAP) )
		f.write( 'FULL_DENS_MAP: %s\n' %(args.full_dens_map) )
		f.write( 'FULL_DENS_MAP_RESO: %0.3f\n' %(args.full_dens_map_reso) )
		f.write( 'SEQ: %s\n' %(SEQ) )
		f.write( 'SS: %s\n' %(SS) )
		f.write( 'NRES: %d\n' %(NRES) )
		f.write( 'NSTRUCT_PER_JOB: %d\n' %(args.nstruct_per_job) )
		f.write( 'CYCLES: %d\n' %(args.cycles) )
		f.write( 'NINIT_POINTS: %d\n' %(NINIT_POINTS) )
		f.write( 'OUT_PREF: %s\n' %(args.out_pref) )
		f.write( 'DISTANCE_CONNECT: %0.3f\n' %(DISTANCE_CONNECT) )
		f.write( 'SHIFT_CENTER: %s\n' %(args.shift_center) )
		f.write( 'PUSH_OUT_FURTHER: %s\n' %(args.push_out_further) )
		f.write( 'NOT_BLIND: %s\n' %(args.not_blind) )
		f.write( 'REPEATS: %d\n' %(args.repeats) )
		f.write( 'MOVE_IN: %s\n' %(args.move_in) )
		f.write( 'USE_ONLY_ONE_END_NODE: %s\n' %(USE_ONLY_ONE_END_NODE) )
		f.write( 'END_NODE_ANGLE_CUTOFF: %0.3f\n' %(END_NODE_ANGLE_CUTOFF) )
		f.write( 'USE_END_NODE: %d\n' %(args.use_end_node) )
		f.write( 'EXTRA_FLAGS: %s\n' %(args.extra_flags) )
		f.write( 'NATIVE: %s\n' %(args.native) )

	
	###################################
	## Step 1: get points from map
	###################################
	
	print( "Converting density map to graph." )
	command_template = 'e2segment3d.py {DENS_MAP} --pdbout={OUT_PREF}_init_points.pdb '
	command_template += '--process=segment.distance:maxsegsep=18:minsegsep=15:thr={MAP_THR}'
	if args.shift_center:
		command_template += ' --shifttocenter'
		formatted_command_template = command_template.format(DENS_MAP=DENS_MAP, OUT_PREF=args.out_pref, MAP_THR=args.map_thr)
		if not args.test:
			os.system( formatted_command_template )
	else:
		formatted_command_template = command_template.format(DENS_MAP=DENS_MAP, OUT_PREF=args.out_pref, MAP_THR=args.map_thr)
		if not args.test:
			os.system( formatted_command_template )
		#print( formatted_command_template )

	if not os.path.exists( '{OUT_PREF}_init_points.pdb'.format( OUT_PREF=args.out_pref ) ):
		print( "ERROR: couldn't convert density map to a graph, please double check your map_thr, it might be too high." )
		exit()
	
	points = []
	pdb_lines = []
	for line in open('{OUT_PREF}_init_points.pdb'.format(OUT_PREF=args.out_pref)):
		if not line.startswith('ATOM'): continue
		x = float(line[29:38].replace(' ',''))
		y = float(line[38:46].replace(' ',''))
		z = float(line[46:54].replace(' ',''))
		pdb_lines.append( line )
		points.append( [x,y,z] )

	# this file is actually useful, let's keep it
	#if not args.debug and not args.test:
	#	os.remove( '{OUT_PREF}_init_points.pdb'.format(OUT_PREF=args.out_pref))
	
	###################################
	## Step 2: Make a graph from the map
	###################################
	G = nx.Graph()
	
	all_distances_btwn = {}
	all_dens_btwn = []
	for i in range(len(points)):
		# add the node
		G.add_node( i, position=np.array(points[i]) )
		for j in range(i+1, len(points)):
			G.add_node( j, position=points[j] )
			# check whether we want to add an edge here
			distance_btwn = get_distance( points[i], points[j] )
			dens_btwn = get_density_btwn_points( DENS_MAP, points[i], points[j], 5 )
			if i not in all_distances_btwn.keys():
				all_distances_btwn[ i ] = {}
			if j not in all_distances_btwn.keys():
				all_distances_btwn[ j ] = {}
			all_distances_btwn[ i ][ j ] = distance_btwn
			all_distances_btwn[ j ][ i ] = distance_btwn
			if distance_btwn < DISTANCE_CONNECT and dens_btwn > 0.02:
				G.add_edge( i, j, weight=get_distance( points[i], points[j] ) )
				all_dens_btwn.append( dens_btwn )
	
	# quickly figure out if we need to add any more connections to "end nodes"
	avg_all_dens_btwn = np.mean( all_dens_btwn )
	for node in G.nodes():
		num_neighbors = get_num_neighbors( G, node )
		if num_neighbors == 1:
			# check if there are any other points within a slightly larger distance
			# but if the distance is larger, then also check that there's lots of density btwn the points
			for node2 in all_distances_btwn[ node ].keys():
				if all_distances_btwn[ node ][ node2 ] < 27.0:
					dens_btwn = get_density_btwn_points( DENS_MAP, G.nodes[node]['position'], G.nodes[node2]['position'], 5 )
					if dens_btwn > avg_all_dens_btwn* 0.5:
						G.add_edge( node, node2, weight=all_distances_btwn[ node ][ node2 ] )
			
	number_of_map_nodes = len( G.nodes() )
	end_nodes_in_map = []
	for node in G.nodes():
		neighbors = get_neighbors( G, node )
		num_neighbors = len( neighbors )
		if num_neighbors == 1:
			end_nodes_in_map.append( node )
		# possibly allow end nodes that have more connections, but have narrow angle between connections
		elif num_neighbors >= 2:
			# maybe this is also an end node
			# check all pairs of 2 and see if the angles are all < cutoff
			possible_end_node = True
			for pair in itertools.combinations( neighbors, 2 ):
				v1 = G.nodes[ pair[0] ]['position'] - G.nodes[node]['position']
				v2 = G.nodes[ pair[1] ]['position'] - G.nodes[node]['position']
				angle = get_angle_btwn( v1, v2 )
				if angle > END_NODE_ANGLE_CUTOFF:
					possible_end_node = False
					break
			if possible_end_node:
				end_nodes_in_map.append( node )
	
	# if there are more than 2 end nodes in the map, then just randomly select one
	# if there are exactly 2, then use both, unless user specified USE_ONLY_ONE_END_NODE
	# or let the use specify the exact end node to use
	print( "Possible end nodes in the map: " )
	end_nodes_str = ''
	for e in end_nodes_in_map:
		end_nodes_str += str(e) + ' '
	print( end_nodes_str )
	print( "You can visualize the end nodes in %s_init_points.pdb" %(args.out_pref ) )
	print( "You can specify which of these end nodes you'd like to use with -use_end_node" )
	if len( end_nodes_in_map ) == 0:
		print( "ERROR: Did not find any end nodes in the map!")
		print( "This can often be fixed by playing with the -map_thr value.")
		exit()
	elif args.use_end_node != -1:
		if args.use_end_node not in end_nodes_in_map:
			print( "The end node you specified is not one of the end nodes found!")
			exit()
		end_nodes_in_map = [ args.use_end_node ]
	elif USE_ONLY_ONE_END_NODE or len( end_nodes_in_map ) > 2:
		end_nodes_in_map = [ end_nodes_in_map[0] ]
	#print( "End nodes: " )
	#print( end_nodes_in_map )
	
	###################################
	## Step 3: Make a graph from the ss
	###################################
	print( "Converting secondary structure to graph." )
	G_ss, stem_connections = DRRAFTER_util.graph_from_ss( SS )
	end_nodes_in_ss = []
	last_res = len( SS ) -1
	first_res = 0
	for node in G_ss.nodes():
		num_neighbors = get_num_neighbors( G_ss, node )
		if num_neighbors == 1:
			end_nodes_in_ss.append( node )
		# also check if this is the first helix and there are no dangling ends
		# that should also be a possible end node
		elif node == 'S0' and first_res in G_ss.nodes[ node ]['residues'] and last_res in G_ss.nodes[ node ]['residues']:
			end_nodes_in_ss.append( node )
	if len( end_nodes_in_ss ) == 0:
		print( "ERROR: Did not find any end nodes in the secondary structure -- auto helix setup will not work for this case" )
		exit()

	if not USE_ONLY_ONE_END_NODE and len( end_nodes_in_ss ) == 1:
		print( "ERROR: Trying to map one secondary structure end node to two end nodes in the map!" )
		print( "The best way to fix this is to play around with the -map_thr until just a single end node is identified in the map." )
		exit()
	
	# need a representative residue number from each element
	# make another graph for the secondary structure to figure out the distance between the elts
	G_res = DRRAFTER_util.residue_graph_from_ss( SS )

	###################################
	## Step 4: Generate possible mappings
	###################################
	print( "Mapping secondary structure to density map." )
	all_possible_mappings = []
	# all end nodes in the map must be filled by end nodes in the SS
	for assignment in itertools.permutations( end_nodes_in_ss, len(end_nodes_in_map) ):
		remaining_end_nodes_ss = [x for x in end_nodes_in_ss if x not in assignment ]
		ss_elt_to_map_elt = {}
		for i, elt in enumerate(assignment):
			ss_elt_to_map_elt[ elt ] = end_nodes_in_map[i]
		all_possible_mappings.append( ss_elt_to_map_elt )
	
	#print( len( all_possible_mappings ), "total mappings" )
	
	
	###################################
	## Step 5: Build RNA helices
	###################################
	print( "Setting up DRRAFTER runs." )
	## also prebuild hairpin loops to save time during the runs
	for node in G_ss.nodes():
		if not node.startswith('S'): continue
		# make the helix
		# if it's a hairpin, then make the hairpin loop
		residues =  G_ss.nodes[ node ][ 'residues' ]
		residues.sort()
		helix_seq = ''
		for rindex, r in enumerate( residues ):
			if rindex > 0:
				if r > residues[ rindex-1 ] +1:
					helix_seq += ' '
					loop_start = residues[ rindex-1] +1
					loop_stop = residues[rindex] -1
			helix_seq += SEQ[r]
	
		tag = DRRAFTER_util.make_tag_with_dashes( residues, ['A' for x in range(len(residues))] )
		elt_name = node.replace('S', 'H')
		helix_pdb_name = '{out_pref}_{name}.pdb'.format( out_pref=args.out_pref, name=elt_name )
		if not os.path.exists( helix_pdb_name ):
			DRRAFTER_util.rna_helix( output_pdb=helix_pdb_name, seq=helix_seq, resnum=tag, 
				rosetta_directory=args.rosetta_directory,rosetta_extension=args.rosetta_extension,
				test=args.test )
		if SS[ loop_start: loop_stop+1].count( '.' ) == loop_stop - loop_start +1:
			if os.path.exists( '{out_pref}_{name}_full.out.1.pdb'.format(out_pref=args.out_pref, name=elt_name)): continue
			# it's a hairpin, so let's prebuild it
			with open( 'fasta_{out_pref}_{name}.txt'.format( out_pref=args.out_pref, name=elt_name), 'w' ) as f:
				f.write( '>{name}.pdb A:{res_first}-{res_last}\n'.format(name=elt_name, res_first=residues[0], res_last=residues[-1]))
				f.write( '{sequence}\n'.format( sequence =SEQ[residues[0]:residues[-1]+1]) )
			with open( 'secstruct_{out_pref}_{name}.txt'.format( out_pref=args.out_pref, name=elt_name), 'w' ) as f:
				f.write( '{ss}\n'.format( ss=SS[residues[0]:residues[-1]+1] ) )
			# flags file
			with open( 'flags_{out_pref}_{name}'.format(out_pref=args.out_pref, name=elt_name), 'w' ) as f:
				f.write( '-fasta fasta_{out_pref}_{name}.txt\n'.format(out_pref=args.out_pref, name=elt_name) )
				f.write( '-secstruct_file secstruct_{out_pref}_{name}.txt\n'.format(out_pref=args.out_pref, name=elt_name) )
				f.write( '-s {helix}\n'.format(helix=helix_pdb_name))
				f.write( '-new_fold_tree_initializer true\n' )
				f.write( '-bps_moves false\n' )
				f.write( '-minimize_rna true\n' )
				f.write( '-out:file:silent {out_pref}_{name}_full.out\n'.format(out_pref=args.out_pref, name=elt_name) )
				f.write( '-cycles 2000\n' )
				f.write( '-nstruct 3\n' )
				f.write( '-ignore_zero_occupancy false\n' )
				f.write( '-use_legacy_job_distributor true\n' )
				f.write( '-no_filters\n' )
				f.write( '-mute core protocols basic\n' )
				if args.test:
					f.write( '-testing:INTEGRATION_TEST\n' )
			command_full_helix = '{rosetta_dir}/rna_denovo{ext} @flags_{out_pref}_{name} > tmp_quiet'.format( rosetta_dir=args.rosetta_directory, 
				ext=args.rosetta_extension, out_pref=args.out_pref, name=elt_name )
			print( "Making full helix %s" %(elt_name) )
			os.system( command_full_helix )
			DRRAFTER_util.extract_lowscore_decoys( '{out_pref}_{name}_full.out'.format(out_pref=args.out_pref,name=elt_name), 1, args.rosetta_directory, args.rosetta_extension, test=args.test )
			#os.system( 'extract_lowscore_decoys.py {out_pref}_{name}_full.out 1 > tmp_quiet'.format(out_pref=args.out_pref,name=elt_name) )
			# remove the files that we don't need, unless in debug mode
			if not args.debug:
				os.remove( 'flags_{out_pref}_{name}'.format(out_pref=args.out_pref, name=elt_name) )
				os.remove( '{out_pref}_{name}_full.out'.format(out_pref=args.out_pref, name=elt_name) )
				os.remove( 'fasta_{out_pref}_{name}.txt'.format( out_pref=args.out_pref, name=elt_name) )
				os.remove( 'secstruct_{out_pref}_{name}.txt'.format( out_pref=args.out_pref, name=elt_name) )
			os.remove( 'tmp_quiet' )
		
	
	###############
	
	flags_template = FLAGS_TEMPLATE
	if args.not_blind:
		flags_template += '-fragment_homology_rmsd 1.2\n-exclusion_match_type MATCH_YR\n'
		flags_template += '-exclude_fragment_files {native}\n'.format(native=args.native)
	
	
	flags_template += args.extra_flags
	
	FASTA_FILE = 'fasta_{out_pref}.txt'.format(out_pref=args.out_pref)
	with open( FASTA_FILE, 'w' ) as f:
		f.write( '>{out_pref}.pdb A:0-{end_num}\n'.format(out_pref=args.out_pref,end_num=len(SEQ)-1))
		f.write( SEQ + '\n' )
	
	SECSTRUCT_FILE = 'secstruct_{out_pref}.txt'.format(out_pref=args.out_pref)
	with open( SECSTRUCT_FILE, 'w') as f:
		f.write( SS + '\n' )
	
	# fit probe helices at each of the end node locations
	# (do this so that we only need to optimize fits once per location,
	# and because the fitting doesn't work very well with short helices)
	
	# make the probe helix if it doesn't exist
	if not os.path.exists( 'probe_helix_6.pdb' ):
		DRRAFTER_util.rna_helix( output_pdb='probe_helix_6.pdb', seq='gggggg cccccc', resnum='A:1-12', 
			rosetta_directory=args.rosetta_directory,rosetta_extension=args.rosetta_extension, test=args.test )
	
	# then figure out all the end nodes in the map that we need to look at
	nodes_for_probe_helices = []
	for mapping in all_possible_mappings:
		# mapping goes from ss_elt : map_elt
		for val in mapping.values(): # going through all the map elts
			nodes_for_probe_helices.append( val )
	
	probe_helices_for_end_nodes_in_map = {}
	for end_node_in_map in set(nodes_for_probe_helices):
		probe_helix_name = 'probe_helix_{num}_{out_pref}.pdb'.format(num=end_node_in_map,out_pref=args.out_pref)
		final_probe_helix_name = 'after_opt_probe_helix_{num}_{out_pref}.pdb'.format(num=end_node_in_map, out_pref=args.out_pref)
		position_of_end_node_in_map = G.nodes[ end_node_in_map ]['position']
		neighbors_map = get_neighbors( G, end_node_in_map )
		position_2 = G.nodes[ neighbors_map[0] ][ 'position' ]
		# push position of end node in map out a little
		if args.push_out_further:
			position_of_end_node_in_map += 2.0*( position_of_end_node_in_map - position_2 )
		elif not args.move_in:
			position_of_end_node_in_map += 0.5*( position_of_end_node_in_map - position_2 )
		# run c++ app to get the aligned helix
		command_align_helix = '{rosetta_directory}/fit_helix_in_map{rosetta_extension} '
		command_align_helix += '-mapfile {map} -end_node_in_map {end_x} {end_y} {end_z} '
		command_align_helix += '-next_node_in_map {next_x} {next_y} {next_z} '
		command_align_helix += '-mute all -probe_helix_pdb probe_helix_6.pdb ' 
		#command_align_helix += '-mute core protocols.moves.RigidBodyMover -probe_helix_pdb probe_helix_6.pdb ' 
		command_align_helix += '-output_pdb_name {out_pdb} -mapreso {map_res} -fit_helix true '
		command_align_helix += '-nrepeats {nrepeats} '
		if args.test:
			command_align_helix += '-testing:INTEGRATION_TEST '
		formatted_command_align_helix = command_align_helix.format( rosetta_directory = args.rosetta_directory, 
			rosetta_extension=args.rosetta_extension,
			map=args.full_dens_map, end_x = position_of_end_node_in_map[0], end_y = position_of_end_node_in_map[1], 
			end_z = position_of_end_node_in_map[2], next_x = position_2[0], next_y = position_2[1],
			next_z = position_2[2], out_pdb=final_probe_helix_name, map_res=args.full_dens_map_reso, nrepeats=args.repeats )
		os.system( formatted_command_align_helix )
		if not args.debug:
			os.remove( 'pose_align_points.pdb' )

		probe_helices_for_end_nodes_in_map[ end_node_in_map ] = final_probe_helix_name
	
	
	# now go through the mappings and do the alignment to the probe helices
	for mapping_index, try_mapping in enumerate(all_possible_mappings):
		ss_elts_from_mapping = try_mapping.keys()
		fit_pdbs = []
		fit_pdbs_not_full = []
		aligned_stems = []
		for ss_elt in ss_elts_from_mapping:
			neighbors = get_neighbors( G_ss, ss_elt )
			if ss_elt.startswith( 'S' ):
				aligned_stems.append( ss_elt )
				residues = G_ss.nodes[ ss_elt ][ 'residues' ]
			else:
				aligned_stems.append( neighbors[0] )
				residues = G_ss.nodes[ neighbors[0] ][ 'residues' ]
			residues.sort()
			helix_seq = ''
			for rindex, r in enumerate(residues):
				if rindex > 0:
					if r > residues[rindex-1] +1:
						helix_seq += ' '
				helix_seq += SEQ[r]
			if ss_elt.startswith( 'H' ):
				helix_pdb_name_1 = '%s_%s.pdb' %(args.out_pref, ss_elt )
			elif ss_elt.startswith( 'S' ):
				helix_pdb_name_1 = '%s_%s.pdb' %(args.out_pref, ss_elt.replace('S','H' ) )
			else:
				helix_pdb_name_1 = '%s_%s.pdb' %(args.out_pref, neighbors[0].replace('S','H') )
			full_helix_pdb_name = ''
			if os.path.exists( '%s_%s_full.out.1.pdb' %(args.out_pref, ss_elt)):
				full_helix_pdb_name = '%s_%s_full.out.1.pdb' %(args.out_pref, ss_elt)
			helices_to_align = [ helix_pdb_name_1 ]
			if full_helix_pdb_name != '': helices_to_align.append( full_helix_pdb_name )
			# do the alignment
			for helix_pdb_name in helices_to_align:
				fit_reverse = False
				# if the end node element is a loop (i.e. the dangling 5' and 3' end of the RNA, then we need to fit the helix in reverse
				if ss_elt.startswith('L') or ss_elt.startswith('S'):
					fit_reverse = True
				# what points are we trying to align to?
				end_node_in_map = try_mapping[ ss_elt ]
				# get the probe helix for this end node
				probe_helix_name = probe_helices_for_end_nodes_in_map[ end_node_in_map ]
				output_aligned_helix_name = 'aligned_map_%d_%s' %( mapping_index, helix_pdb_name )
	
				### do alignment in Rosetta
				command_alignment = '{rosetta_dir}/fit_helix_in_map{ext} -mute all -align true '
				#command_alignment = '{rosetta_dir}/fit_helix_in_map{ext} -mute core -align true '
				command_alignment += '-helix_for_alignment {helix_for_alignment} -fit_probe_helix {fit_probe_helix} '
				command_alignment += '-output_pdb_name_aligned {output_pdb} '
				if args.test:
					command_alignment += '-testing:INTEGRATION_TEST '
				formatted_command_alignment = command_alignment.format( rosetta_dir=args.rosetta_directory, ext=args.rosetta_extension, 
						helix_for_alignment = helix_pdb_name, fit_probe_helix=probe_helix_name, 
						output_pdb=output_aligned_helix_name )
				if fit_reverse:
					formatted_command_alignment += '-fit_reverse true'
					#print( formatted_command_alignment )
					os.system( formatted_command_alignment )
				else:
					command_alignment += '-fit_reverse false'
					#print( formatted_command_alignment )
					os.system( formatted_command_alignment )

			if full_helix_pdb_name != '':
				fit_pdbs.append( 'aligned_map_%d_%s' %( mapping_index, full_helix_pdb_name ) )
			else:
				fit_pdbs.append( 'aligned_map_%d_%s' %( mapping_index, helix_pdb_name_1 ) )
			fit_pdbs_not_full.append( 'aligned_map_%d_%s' %( mapping_index, helix_pdb_name_1 ) )
	
		###################################
		## Step 5: Set up DRRAFTER run
		###################################
		fit_pdbs_str = ''
		for pdb in fit_pdbs:
			fit_pdbs_str += pdb
			fit_pdbs_str += ' '
		cat_all_aligned_pdb = 'all_aligned_%s_%d.pdb' %( args.out_pref, mapping_index )
		os.system( 'cat %s > %s' %( fit_pdbs_str, cat_all_aligned_pdb ))
		DRRAFTER_util.better_reorder_pdb( cat_all_aligned_pdb )
		cat_all_aligned_pdb_reorder = 'all_aligned_%s_%d.REORDER.pdb' %(args.out_pref, mapping_index )
		# delete the old files (before reordering)
		os.remove( cat_all_aligned_pdb )
		# the other helices...
		other_helices = ''
		other_helices_no_full = ''
		for node in G_ss.nodes():
			if node in aligned_stems: continue
			if not node.startswith( 'S' ): continue
			helix_pdb_name = '%s_%s.pdb' %(args.out_pref, node.replace('S','H') )
			full_pdb_name = '%s_%s_full.out.1.pdb' %(args.out_pref, node.replace('S','H') )
			other_helices_no_full += helix_pdb_name
			other_helices_no_full += ' '
			if os.path.exists( full_pdb_name ):
				other_helices += full_pdb_name
				other_helices += ' '
			else:
				other_helices += helix_pdb_name
				other_helices += ' '
	
		with open( 'flags_%s_%d_R1' %(args.out_pref, mapping_index), 'w') as f:
			text = flags_template.format( 
				cat_helix = cat_all_aligned_pdb_reorder, other_helices = other_helices, 
				out_pref = args.out_pref + '_' + str(mapping_index), 
				fasta = FASTA_FILE, secstruct = SECSTRUCT_FILE, map_file = args.full_dens_map, 
				map_reso=args.full_dens_map_reso, nstruct=args.nstruct_per_job, cycles=args.cycles )
			f.write( text )
			if args.test:
				f.write( '-testing:INTEGRATION_TEST\n' )
	
		with open( 'command_%s_%d_R1' %(args.out_pref, mapping_index ), 'w' ) as f:
			f.write( '%s/rna_denovo%s @flags_%s_%d_R1\n' %( args.rosetta_directory, args.rosetta_extension, args.out_pref, mapping_index ) )
		if not args.debug:
			for pdb in fit_pdbs:
				if os.path.exists( pdb ): os.remove( pdb )
			for pdb in fit_pdbs_not_full:
				if os.path.exists( pdb ): os.remove( pdb )

	# write a file with all of the fits
	with open( '%s_auto_fits.txt' %(args.out_pref ), 'w' ) as f:
		for mapping_number in range(len(all_possible_mappings)):
			f.write( '%s\n' %(mapping_number) )
	
	# clean up
	if not args.debug:
		os.remove( 'probe_helix_6.pdb' )
		if os.path.exists( 'before_opt.pdb' ):
			os.remove( 'before_opt.pdb' )
		for x in probe_helices_for_end_nodes_in_map.keys():
			os.remove( probe_helices_for_end_nodes_in_map[x] )
		if not args.test:
			os.remove( DENS_MAP )

if __name__ == '__main__':
	parser = argparse.ArgumentParser( description='Auto setup helices for DRRAFTER-RNA' )
	parser.add_argument( '-map_thr', type=float, required=True )
	parser.add_argument( '-full_dens_map', type=str, required=True )
	parser.add_argument( '-full_dens_map_reso', type=float, required=True )
	parser.add_argument( '-fasta', type=str, required=True )
	parser.add_argument( '-secstruct', type=str, required=True )
	parser.add_argument( '-out_pref', type=str, required=True )
	parser.add_argument( '-shift_center', action='store_true', default=False )
	parser.add_argument( '-push_out_further', action='store_true', default=False )
	parser.add_argument( '-not_blind', action='store_true', default=False )
	parser.add_argument( '-debug', action='store_true', default=False )
	parser.add_argument( '-test', action='store_true', default=False, help="for running tests on the testing server" )
	parser.add_argument( '-move_in', action='store_true', default=False )
	parser.add_argument( '-repeats', type=int, default=10 )
	parser.add_argument( '-nstruct_per_job', type=int, default=200 )
	parser.add_argument( '-cycles', type=int, default=30000 )
	parser.add_argument( '-use_end_node', type=int, default=-1 )
	parser.add_argument( '-extra_flags', type=str, default='' )
	parser.add_argument( '-native', type=str, default='' )
	parser.add_argument( '-just_low_pass', action='store_true', default=False )
	parser.add_argument( '-fit_only_one_helix', action='store_true', default=False )
	parser.add_argument( '-rosetta_directory', required=True, help="Path to Rosetta executables" )
	parser.add_argument( '-rosetta_extension', default='', help="Extension for Rosetta executables e.g. '.linuxgccrelease'" )

	args = parser.parse_args()
	auto_setup_helices( args )
