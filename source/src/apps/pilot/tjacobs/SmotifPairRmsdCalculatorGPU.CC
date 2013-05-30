// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

#include <devel/sewing/BundlePairRmsdWorkUnit.hh>

#include <basic/gpu/GPU.hh>
#include <basic/gpu/Timer.hh>

#include <core/pose/util.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

#include <cstdio>

#include <devel/init.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

// ObjexxFCL libraries
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//TEST
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("BundlePairRmsdCalculator");

static int BLOCK_SIZE = 512;

namespace BundlePairRmsdCalculator {
	basic::options::IntegerOptionKey const total_fractions( "total_fractions" ); // fraction to run
	basic::options::IntegerOptionKey const fraction( "fraction" ); // fraction to run
	basic::options::IntegerOptionKey const subfraction_size( "subfraction_size" ); // subfraction size
}

struct gpu_node_block_ptrs
{
	cl_mem d_bundle_ids; // pointer to an array of integers representing the bundle ids for the set of nodes in this block
	cl_mem coords_for_rms_calc; // pointer to an array of xyzvectors that should be used to calculate RMDS
	cl_mem rms_coords_moment_of_inertia; // pointer to an array of floats, one per node, representing the moment of inertia of the atoms in the rms calculation after they have been translated so that their center of mass is at the origin.
	cl_mem coords_for_col_calc; // pointer to an array of xyzvectors for clash-checking purposes
	cl_mem coords_for_col_calc_offset; // pointer to an array of offsets representing the first coordinate in the coords_for_col_calc belonging to a given node
	cl_mem n_coords_for_col_calc; // pointer to an array of integers representing the number of coordinates in belonging to each node in the coords_for_col_cal array
	cl_mem comparison_coord_ind_list; // pointer to an array integers representing those atoms for each node which should be clash checked against
	cl_mem comparison_coord_ind_offset; // pointer to an array of offsets representing the beginning index into the comparison_coord_ind_list for each node
	cl_mem n_comparison_coords; // pointer to an array of integers representing the number of comparison coordinates for each node
};

struct node_data
{
	int node_index;
	int bundle_index;
	int helix_1_ind;
	int helix_2_ind;
	int n_atoms_in_rms_calc;
};

struct calculation_data_for_one_block_pair
{
	ObjexxFCL::FArray2D< float > rms_data_table;
	ObjexxFCL::FArray2D< float > clash_data_table;
	utility::vector1< float > rms_data;
	utility::vector1< float > clash_data;
	utility::vector1< int > good_node_pairs;
	int n_good_node_pairs;
	cl_mem d_rmsd_table;
	cl_mem d_clash_table;
	cl_mem d_calculate_pair_table;
	cl_mem d_good_node_pairs;
	cl_mem d_n_good_node_pairs;
	cl_event read_n_good_node_pairs_event;
};

typedef std::map<core::Size, utility::vector1< numeric::xyzVector<core::Real> > > helix_coord_map;
typedef std::map<core::Size, utility::vector1< core::Size > > bundle_helices_map;
typedef std::map<core::Size, bool > helix_flipped_map;

std::vector< std::string >
gpu_rmsd_programs() {
	std::vector< std::string > programs;
	programs.push_back( "/home/tjacobs2/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_xyzfunctions.cl" );
	programs.push_back( "/home/tjacobs2/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_rmsd_functions.cl" );
	programs.push_back( "/home/tjacobs2/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_helical_bundle_rms_and_clash_calculations.cl" );
	return programs;
}

void
load_coordinates_from_database(
	utility::sql_database::sessionOP db_session,
	helix_coord_map & helix_coords,
	bundle_helices_map & bundle_helices,
	helix_flipped_map & helix_flipped
)
{
	/**Get coordinates**/
	std::string select_helix_coords =
//	"SELECT bh.bundle_id, bh.helix_id, bh.flipped, res.x, res.y, res.z\n"
//	"FROM bundle_helices bh\n"
//	"JOIN structures s ON\n"
//	"	s.struct_id = bh.struct_id\n"
//	"JOIN residue_atom_coords res ON\n"
//	"	bh.struct_id = res.struct_id AND\n"
//	"	res.atomno IN (1,2,3,4) AND\n"
//	"	(res.seqpos BETWEEN bh.residue_begin AND bh.residue_end)\n"
//	"ORDER BY bh.helix_id, res.seqpos, res.atomno;";
	
	"SELECT sp.smotif_pair_id, sss.segment_id, res.x, res.y, res.z\n"
	"FROM smotif_pairs sp ON\n"
	"JOIN smotifs s1 ON\n"
	"	sp.smotif_id_1 = s1.smotif_id\n"
	"JOIN smotifs s2 ON\n"
	"	sp.smotif_id_2 = s2.smotif_id\n"
	"JOIN secondary_structure_segments sss ON\n"
	"	s1.secondary_struct_segment_id_1 = sss.segment_id OR\n"
	"	s1.secondary_struct_segment_id_2 = sss.segment_id OR\n"
	"	s2.secondary_struct_segment_id_1 = sss.segment_id OR\n"
	"	s2.secondary_struct_segment_id_2 = sss.segment_id\n"
	"JOIN residue_atom_coords res ON\n"
	"	sm.struct_id = res.struct_id AND\n"
	"	res.atomno IN (1,2,3,4) AND\n"
	"	(res.seqpos BETWEEN sss.residue_begin AND sss.residue_end)\n"
	"ORDER BY sss.segment_id, res.seqpos, res.atomno;";
	
	cppdb::statement coords_stmt=basic::database::safely_prepare_statement(select_helix_coords, db_session);
	cppdb::result coords_res=basic::database::safely_read_from_database(coords_stmt);
	TR << "Done querying coordinates from DB" << std::endl;

	while(coords_res.next())
	{
		core::Size bundle_id, helix_id;
		core::Real x, y, z;

		coords_res >> bundle_id >> helix_id >> x >> y >> z;
		//std::cout << "bundle id: " << bundle_id << " helix_id " << helix_id << " " << x << " " << y << " " << z << "\n";
		helix_coords[helix_id].push_back(numeric::xyzVector< core::Real >(x,y,z));
		if ( std::find( bundle_helices[ bundle_id ].begin(), bundle_helices[ bundle_id ].end(), helix_id ) == bundle_helices[ bundle_id ].end() ) bundle_helices[bundle_id].push_back(helix_id);
		//helix_flipped[helix_id]=flipped;
	}
	//std::cout << std::endl;
}

void
ship_coordinates_to_gpu(
	basic::gpu::GPU & gpu,
	helix_coord_map const & helix_coords,
	bundle_helices_map const & bundle_helices,
	utility::vector1< node_data > & ndat,
	core::Size blockid,
	gpu_node_block_ptrs & block_ptrs
)
{
	// ok -- count up how large each comparison helix is, and how large each clash-check helix is
	// then allocate a set of arrays on the CPU for their coordinates and their offsets,
	// then allocate space on the GPU for these arrays and send the coordinates to the GPU.
	core::Size n_atoms_in_rms_helices = 0;
	core::Size count_atoms_in_col_helices = 0;
	core::Size count_n_col_helices = 0;

	utility::vector1< int > n_col_coords( BLOCK_SIZE, 0 );
	utility::vector1< int > col_coords_offsets( BLOCK_SIZE, 0 );

	utility::vector1< int > n_comparison_coords( BLOCK_SIZE, 0 );
	utility::vector1< int > comparison_coords_offsets( BLOCK_SIZE, 0 );

	int num_nodes_this_block = (blockid)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (blockid-1)*BLOCK_SIZE;

	utility::vector1< int > block_ids_for_nodes( BLOCK_SIZE, 0 );
	for ( core::Size ii = 1; ii <= (core::Size) num_nodes_this_block; ++ii ) {
		block_ids_for_nodes[ ii ] = ndat[ii+(blockid-1)*BLOCK_SIZE].bundle_index;
	}

	for ( core::Size ii = 1; ii <= (core::Size) num_nodes_this_block; ++ii ) {
		core::Size const iinode = ii + (blockid-1)*BLOCK_SIZE;
		assert( iinode <= ndat.size() );

		core::Size ii_bundle = ndat[iinode].bundle_index;
		utility::vector1< core::Size > const & iihelices = bundle_helices.find( ii_bundle )->second;
		//if ( iinode <= 5 ) {
		//	for ( core::Size jj = 1; jj <= iihelices.size(); ++jj ) {
		//		std::cout << "node " << ii << " helix " << iihelices[jj] << std::endl;
		//	}
		//}

		core::Size ii_n_atoms_in_rms_helices = 0;
		core::Size ii_n_atoms_in_col_helices = 0;
		core::Size ii_n_col_helices = 0;
		for ( core::Size jj = 1; jj <= iihelices.size(); ++jj ) {
			core::Size jjhelix_index = iihelices[jj];
			utility::vector1< numeric::xyzVector< core::Real> > const & jjhelix = helix_coords.find( jjhelix_index )->second;
			assert( jjhelix.size() != 0 );
			if ( (int) jjhelix_index == ndat[iinode].helix_1_ind || (int) jjhelix_index == ndat[iinode].helix_2_ind ) {
				ii_n_atoms_in_rms_helices += jjhelix.size();
			} else {
				++ii_n_col_helices;
				ii_n_atoms_in_col_helices += jjhelix.size();
			}
		}

		if ( n_atoms_in_rms_helices == 0 ) {
			n_atoms_in_rms_helices = ii_n_atoms_in_rms_helices;
		} else if ( ii_n_atoms_in_rms_helices != n_atoms_in_rms_helices ) {
			// error!
			utility_exit_with_message( "Sanity check failed. Node " + utility::to_string( ii ) + " with RMS helices of index " +
				utility::to_string( ndat[iinode].helix_1_ind ) + " and " + utility::to_string( ndat[iinode].helix_2_ind ) + " have more atoms " +
				"in them than the other helices so far encountered; " + utility::to_string( n_atoms_in_rms_helices ) + " vs " +
				utility::to_string( ii_n_atoms_in_rms_helices ));
		}
		ndat[iinode].n_atoms_in_rms_calc = ii_n_atoms_in_rms_helices;
		count_atoms_in_col_helices += ii_n_atoms_in_col_helices;
		count_n_col_helices += ii_n_col_helices;
		n_col_coords[ ii ] = ii_n_atoms_in_col_helices;
		if ( ii > 1 ) col_coords_offsets[ ii ] = col_coords_offsets[ii-1] + n_col_coords[ii-1];
		n_comparison_coords[ ii ] = ii_n_col_helices; // <-- pick one coordinate for each helix to clash check against all atoms in the comparison helices
		if ( ii > 1 ) comparison_coords_offsets[ ii ] = comparison_coords_offsets[ ii-1 ] + n_comparison_coords[ ii-1 ];
		//if ( ii <= 20 ) {
		//	std::cout << "node " << ii << ":\n";
		//	std::cout << " ii_n_atoms_in_rms_helices: " << ii_n_atoms_in_rms_helices << "\n";
		//	std::cout << " ii_n_col_helices: " << ii_n_col_helices << "\n";
		//	std::cout << " ii_n_atoms_in_col_helices: " << ii_n_atoms_in_col_helices << "\n";
		//	std::cout << " col_coords_offsets[ ii ]: " << col_coords_offsets[ii] << "\n";
		//	std::cout << " comparison_coords_offsets[ ii ]: " << comparison_coords_offsets[ ii ] << "\n";
		//}
	}

	/// allocate space for the coordinates
	utility::vector1< numeric::xyzVector< float > > all_rms_coords( n_atoms_in_rms_helices * BLOCK_SIZE );
	utility::vector1< numeric::xyzVector< float > > all_col_coords( count_atoms_in_col_helices );
	utility::vector1< int > all_comparison_coord_inds( count_n_col_helices );

	core::Size count_rms_coords = 0, count_col_coords = 0, count_comp_coords = 0;
	for ( core::Size ii = 1; ii <= (core::Size) num_nodes_this_block; ++ii ) {
		core::Size iinode = ii+(blockid-1)*BLOCK_SIZE;
		assert( iinode <= ndat.size() );

		core::Size ii_bundle = ndat[iinode].bundle_index;
		utility::vector1< core::Size > const & iihelices = bundle_helices.find( ii_bundle )->second;
		core::Size jj_col_helices_atom_offset = 0;
		{ // scope  -- add helix 1, then add helix 2
			utility::vector1< numeric::xyzVector< core::Real> > const & helix1 = helix_coords.find( ndat[iinode].helix_1_ind )->second;
			for ( core::Size kk = 1; kk <= helix1.size(); ++kk ) {
				++count_rms_coords;
				assert( count_rms_coords <= n_atoms_in_rms_helices * BLOCK_SIZE );
				//if ( ii == 181 || ii == 196 ) {
				//	std::cout << "coords: " << ii << " " << helix1[kk].x() << " " << helix1[kk].y() << " " << helix1[kk].z() << std::endl;
				//}
				all_rms_coords[ count_rms_coords ] = helix1[kk];
			}
			utility::vector1< numeric::xyzVector< core::Real> > const & helix2 = helix_coords.find( ndat[iinode].helix_2_ind )->second;
			for ( core::Size kk = 1; kk <= helix2.size(); ++kk ) {
				++count_rms_coords;
				assert( count_rms_coords <= n_atoms_in_rms_helices * BLOCK_SIZE );
				//if ( ii == 181 || ii == 196 ) {
				//	std::cout << "coords: " << ii << " " << helix2[kk].x() << " " << helix2[kk].y() << " " << helix2[kk].z() << std::endl;
				//}
				all_rms_coords[ count_rms_coords ] = helix2[kk];
			}
		}

		for ( core::Size jj = 1; jj <= iihelices.size(); ++jj ) {
			//if ( ii <= 20 ) { std::cout << "iihelices: " << ii << " " << jj << " " << count_rms_coords << std::endl; }
			core::Size jjhelix_index = iihelices[jj];
			utility::vector1< numeric::xyzVector< core::Real> > const & jjhelix = helix_coords.find( jjhelix_index )->second;
      assert( jjhelix.size() != 0 );
      if ( (int) jjhelix_index == ndat[iinode].helix_1_ind || (int) jjhelix_index == ndat[iinode].helix_2_ind ) {
				continue;
      } else {
				for (core::Size kk = 1; kk <= jjhelix.size(); ++kk ) {
					++count_col_coords;
					assert( count_col_coords <= count_atoms_in_col_helices );
					all_col_coords[ count_col_coords ] = jjhelix[kk];
					//if ( ii == 181 || ii == 196 ) {
					//	std::cout << "clashcoords: " << ii << " " << jjhelix[kk].x() << " " << jjhelix[kk].y() << " " << jjhelix[kk].z() << std::endl;
					//}
				}
				assert( jjhelix.size() % 4 == 0 ); // there should be four atoms per residue and whole residues only in this structure.
				++count_comp_coords;
				assert( count_comp_coords <= count_n_col_helices );
				// indicate C-alpha in 0-based indexing on the resude in the center of the helix + offset from previous colision helices
				core::Size jj_comp_coord_index = (( jjhelix.size() / 4 ) / 2) * 4 + 1 + jj_col_helices_atom_offset;
				//if ( ii <= 20 ) { std::cout << " jj_comp_coord_index: " << ii << " " << jj << " " << jj_comp_coord_index << std::endl; }
				all_comparison_coord_inds[ count_comp_coords ] = jj_comp_coord_index;
				jj_col_helices_atom_offset += jjhelix.size();
      }
		}
	}

	cl_mem d_bundle_ids = gpu.AllocateMemory( sizeof( int ) * BLOCK_SIZE, &block_ids_for_nodes[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_bundle_ids ) utility_exit_with_message( "Failed to allocate memory on the GPU for bundle_ids");

	cl_mem d_rms_coords = gpu.AllocateMemory( sizeof( float ) * 3 * n_atoms_in_rms_helices * BLOCK_SIZE, &all_rms_coords[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_rms_coords ) utility_exit_with_message( "Failed to allocate memory on the GPU for d_rms_coords");

	cl_mem d_rms_coords_moment_of_inertia = gpu.AllocateMemory( sizeof( float ) * BLOCK_SIZE, NULL, CL_MEM_READ_WRITE );
	if ( ! d_rms_coords_moment_of_inertia ) utility_exit_with_message( "failed to allocate memory on the GPU for d_rms_coords_moment_of_inertia" );

	cl_mem d_col_coords = gpu.AllocateMemory( sizeof( float ) * 3 * count_atoms_in_col_helices, &all_col_coords[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_col_coords ) utility_exit_with_message( "Failed to allocate memory on the GPI for d_col_coords");

	cl_mem d_col_coords_offsets = gpu.AllocateMemory( sizeof( int ) * BLOCK_SIZE, &col_coords_offsets[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_col_coords_offsets ) utility_exit_with_message( "Failed to allocate memory on the GPI for d_col_coords_offsets");

 	cl_mem d_n_col_coords = gpu.AllocateMemory( sizeof( int ) * BLOCK_SIZE, &n_col_coords[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_n_col_coords ) utility_exit_with_message( "Failed to allocate memory on the GPI for d_n_col_coords");

 	cl_mem d_comparison_coord_inds = gpu.AllocateMemory( sizeof( int ) * count_n_col_helices, &all_comparison_coord_inds[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_comparison_coord_inds ) utility_exit_with_message( "Failed to allocate memory on the GPI for d_comparison_coord_inds");

	cl_mem d_n_comparison_coords = gpu.AllocateMemory( sizeof( int ) * BLOCK_SIZE, &n_comparison_coords[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_n_comparison_coords ) utility_exit_with_message( "Failed to allocate memory on the GPI for d_n_comparison_coords");

	cl_mem d_comparison_coords_offsets = gpu.AllocateMemory( sizeof( int ) * BLOCK_SIZE, &comparison_coords_offsets[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
	if ( ! d_comparison_coords_offsets ) utility_exit_with_message( "Failed to allocate memory on the GPI for d_comparison_coords_offsets");

	block_ptrs.d_bundle_ids                   = d_bundle_ids;
	block_ptrs.coords_for_rms_calc          = d_rms_coords;
	block_ptrs.rms_coords_moment_of_inertia = d_rms_coords_moment_of_inertia;
	block_ptrs.coords_for_col_calc          = d_col_coords;
	block_ptrs.coords_for_col_calc_offset   = d_col_coords_offsets;
	block_ptrs.n_coords_for_col_calc        = d_n_col_coords;
	block_ptrs.comparison_coord_ind_list    = d_comparison_coord_inds;
	block_ptrs.comparison_coord_ind_offset  = d_comparison_coords_offsets;
	block_ptrs.n_comparison_coords          = d_n_comparison_coords;

	// OK: now translate the rms coordinates so that their center of mass is at the origin and
	// then translate the collision coordinates using the same translation vector as the rms coordinates

	//std::cout << "translate_rmscoords_to_origin: " << num_nodes_this_block << " " << n_atoms_in_rms_helices << std::endl;

	gpu.ExecuteKernel( "translate_rmscoords_to_origin", BLOCK_SIZE, 32, 32,
		GPU_INT, num_nodes_this_block,
		GPU_INT, n_atoms_in_rms_helices,
		GPU_DEVMEM, block_ptrs.coords_for_rms_calc,
		GPU_DEVMEM, block_ptrs.rms_coords_moment_of_inertia,
		GPU_DEVMEM, block_ptrs.coords_for_col_calc,
		GPU_DEVMEM, block_ptrs.coords_for_col_calc_offset,
		GPU_DEVMEM, block_ptrs.n_coords_for_col_calc,
		NULL );
}

void
compute_rms_and_clash_score_for_block_pair(
	basic::gpu::GPU & gpu,
	utility::vector1< node_data > const & ndat,
	core::Size const block1,
	core::Size const block2,
	gpu_node_block_ptrs const & block1ptrs,
	gpu_node_block_ptrs const & block2ptrs,
	float maximum_rmsd_tolerance,
	float minimum_square_distance_tolerance,
	calculation_data_for_one_block_pair & bpd // block-pair-data
)
{

	int block1_nnodes = (block1)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (block1-1)*BLOCK_SIZE;
	int block2_nnodes = (block2)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (block2-1)*BLOCK_SIZE;
	int block1_offset = (block1-1)*BLOCK_SIZE;
	//int block2_offset = (block2-1)*BLOCK_SIZE;

	// allocate space on the GPU for these three tables if that space has not already been allocated
	// Note that the space being used for this block1/block2 combo should definitely no longer
	// be in use for some other block1/block2 combo!
	if ( ! bpd.d_rmsd_table ) bpd.d_rmsd_table = gpu.AllocateMemory( sizeof( float ) * BLOCK_SIZE * BLOCK_SIZE, NULL, CL_MEM_READ_WRITE );
	if ( ! bpd.d_rmsd_table ) utility_exit_with_message( "Failed to allocate memory on the GPU for d_rmsd_table");

	if ( ! bpd.d_clash_table ) bpd.d_clash_table = gpu.AllocateMemory( sizeof( float ) * BLOCK_SIZE * BLOCK_SIZE, NULL, CL_MEM_READ_WRITE );
	if ( ! bpd.d_clash_table ) utility_exit_with_message( "Failed to allocate memory on the GPU for d_clash_table");

	if ( ! bpd.d_calculate_pair_table ) bpd.d_calculate_pair_table = gpu.AllocateMemory( sizeof( unsigned char ) * BLOCK_SIZE * BLOCK_SIZE, NULL, CL_MEM_READ_WRITE );
	if ( ! bpd.d_calculate_pair_table ) utility_exit_with_message( "Failed to allocate memory on the GPU for d_calculate_pair_table");

	if ( ! bpd.d_good_node_pairs ) bpd.d_good_node_pairs = gpu.AllocateMemory( sizeof( int ) * BLOCK_SIZE * BLOCK_SIZE, NULL, CL_MEM_READ_WRITE );
	if ( ! bpd.d_good_node_pairs ) utility_exit_with_message( "Failed to allocate memory on the GPU for d_good_node_pairs");

	if ( ! bpd.d_n_good_node_pairs ) bpd.d_n_good_node_pairs = gpu.AllocateMemory( sizeof( int ), NULL, CL_MEM_READ_WRITE );
	if ( ! bpd.d_n_good_node_pairs ) utility_exit_with_message( "Failed to allocate memory on the GPU for d_n_good_node_pairs");

	cl_kernel calc_rmsd_and_clash_score_kernel = gpu.BuildKernel( "compute_rmsd_and_clash_scores" );

	cl_int errNum;
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 0, sizeof(int), & BLOCK_SIZE );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 0, BLOCK_SIZE" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 1, sizeof(int), & block1_nnodes );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 1, block1_nnodes" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 2, sizeof(int), & block2_nnodes );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 2, block2_nnodes" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 3, sizeof(int), & ndat[1+block1_offset].n_atoms_in_rms_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 3, ndat[ii+block1_offset].n_atoms_in_rms_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 4, sizeof(cl_mem), & block1ptrs.d_bundle_ids );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 4, block1ptrs.d_bundle_ids" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 5, sizeof(cl_mem), & block2ptrs.d_bundle_ids );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, block2ptrs.d_bundle_ids" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 6, sizeof(cl_mem), & block1ptrs.coords_for_rms_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 4, block1ptrs.coords_for_rms_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 7, sizeof(cl_mem), & block2ptrs.coords_for_rms_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, block2ptrs.coords_for_rms_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 8, sizeof(cl_mem), & block1ptrs.rms_coords_moment_of_inertia );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 6, block1ptrs.rms_coords_moment_of_inertia" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel, 9, sizeof(cl_mem), & block2ptrs.rms_coords_moment_of_inertia );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 7, block2ptrs.rms_coords_moment_of_inertia" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,10, sizeof(cl_mem), & block1ptrs.coords_for_col_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 8, block1ptrs.coords_for_col_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,11, sizeof(cl_mem), & block2ptrs.coords_for_col_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 9, block2ptrs.coords_for_col_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,12, sizeof(cl_mem), & block1ptrs.coords_for_col_calc_offset );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 10, block1ptrs.coords_for_col_calc_offset" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,13, sizeof(cl_mem), & block2ptrs.coords_for_col_calc_offset );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 11, block2ptrs.coords_for_col_calc_offset" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,14, sizeof(cl_mem), & block1ptrs.n_coords_for_col_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 12, block1ptrs.n_coords_for_col_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,15, sizeof(cl_mem), & block2ptrs.n_coords_for_col_calc );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 13, block2ptrs.n_coords_for_col_calc" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,16, sizeof(cl_mem), & block1ptrs.comparison_coord_ind_list );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 14, block1ptrs.comparison_coord_ind_list" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,17, sizeof(cl_mem), & block2ptrs.comparison_coord_ind_list );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 15, block2ptrs.comparison_coord_ind_list" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,18, sizeof(cl_mem), & block1ptrs.comparison_coord_ind_offset );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 16, block1ptrs.comparison_coord_ind_offset" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,19, sizeof(cl_mem), & block2ptrs.comparison_coord_ind_offset );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 17, block2ptrs.comparison_coord_ind_offset" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,20, sizeof(cl_mem), & block1ptrs.n_comparison_coords );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 18, block1ptrs.n_comparison_coords" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,21, sizeof(cl_mem), & block2ptrs.n_comparison_coords );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 19, block2ptrs.n_comparison_coords" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,22, sizeof(cl_mem), & bpd.d_rmsd_table );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 20, d_rmsd_table" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,23, sizeof(cl_mem), & bpd.d_clash_table );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 21, d_clash_table" ); }
	errNum = clSetKernelArg(calc_rmsd_and_clash_score_kernel,24, sizeof(cl_mem), & bpd.d_calculate_pair_table );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 22, d_calculate_pair_table" ); }

	basic::gpu::GPU_DEV const & gpu_dev = gpu.device();

	size_t global_work_size[2] = { BLOCK_SIZE, BLOCK_SIZE };
	size_t local_work_size[2] = {8, 8};
	cl_event calcRMSDAndClashKernelEvent;

	//std::cout << "calc_rmsd_and_clash_score_kernel: " << block1 << " " << block2 << " " << ndat[1+block1_offset].n_atoms_in_rms_calc << std::endl;

	//basic::gpu::Timer timer;
	errNum = clEnqueueNDRangeKernel( gpu_dev.commandQueue, calc_rmsd_and_clash_score_kernel,
		2, NULL, global_work_size, local_work_size, 0, NULL, &calcRMSDAndClashKernelEvent );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to run calc rms and clash score kernel" ); }

	// DEBUG: pull back all the results and put them into an FArray *BEFORE* running scan, which would overwrite the first section of the array
	//cl_event read_events[2];
	//errNum = clEnqueueReadBuffer( gpu_dev.commandQueue, bpd.d_rmsd_table, CL_FALSE, 0, sizeof( float )*BLOCK_SIZE*BLOCK_SIZE, &bpd.rms_data_table(1,1), 1, &calcRMSDAndClashKernelEvent, &read_events[0] );
	//errNum = clEnqueueReadBuffer( gpu_dev.commandQueue, bpd.d_clash_table, CL_FALSE, 0, sizeof( float )*BLOCK_SIZE*BLOCK_SIZE, &bpd.clash_data_table(1,1), 1, &calcRMSDAndClashKernelEvent, &read_events[1] );

	// now set up a follow-up kernel launch which will sort through the rms/clash scores and select only the good ones
	// so that they can be read back to the CPU
	cl_kernel scan_good_rmsd_clash_pairs_kernel = gpu.BuildKernel( "scan_good_rmsd_clash_pairs" );

  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 0, sizeof(int), & BLOCK_SIZE );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 0, BLOCK_SIZE" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 1, sizeof(int), & block1_nnodes );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 1, block1_nnodes" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 2, sizeof(int), & block2_nnodes );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 2, block2_nnodes" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 3, sizeof(float), & maximum_rmsd_tolerance );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 3, maximum_rmsd_tolerance" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 4, sizeof(float), & minimum_square_distance_tolerance );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 4, block1ptrs.d_bundle_ids" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 5, sizeof(cl_mem), & bpd.d_calculate_pair_table );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, d_calculate_pair_table" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 6, sizeof(cl_mem), & bpd.d_rmsd_table );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, d_rmsd_table" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 7, sizeof(cl_mem), & bpd.d_clash_table );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 4, d_clash_table" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 8, sizeof(cl_mem), & bpd.d_good_node_pairs );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, d_good_node_pairs" ); }
  errNum = clSetKernelArg(scan_good_rmsd_clash_pairs_kernel, 9, sizeof(cl_mem), & bpd.d_n_good_node_pairs );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, d_n_good_node_pairs" ); }

	cl_event scanKernelEvent;
	size_t two_fifty_six = 256;
	errNum = clEnqueueNDRangeKernel( gpu_dev.commandQueue, scan_good_rmsd_clash_pairs_kernel,
    1, NULL, &two_fifty_six, &two_fifty_six, 1, &calcRMSDAndClashKernelEvent, &scanKernelEvent );
	//errNum = clEnqueueNDRangeKernel( gpu_dev.commandQueue, scan_good_rmsd_clash_pairs_kernel,
  //  1, NULL, &two_fifty_six, &two_fifty_six, 2, read_events, &scanKernelEvent );
	if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to enqueue scan_good_rmsd_clash_pairs_kernel" ); }

	// Now enqueue a read event for the single integer representing the number of good indices
  cl_event read_n_good_pairs_event;
  errNum = clEnqueueReadBuffer( gpu_dev.commandQueue, bpd.d_n_good_node_pairs,  CL_FALSE, 0, sizeof( int ), &bpd.n_good_node_pairs, 1, &scanKernelEvent, &read_n_good_pairs_event );
  if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to enqueue ReadBuffer for n_good_node_pairs" ); }

	/// Wait until this event completes and then proceed to extract data from
	bpd.read_n_good_node_pairs_event = read_n_good_pairs_event;

}

void
retrieve_rms_and_clash_scores_from_gpu(
	basic::gpu::GPU & gpu,
	calculation_data_for_one_block_pair & bpd // block-pair data
)
{
	basic::gpu::GPU_DEV const & gpu_dev = gpu.device();

	std::cout << "retrieving data for " << bpd.n_good_node_pairs << " good node pairs" << std::endl;
	cl_event read_buffer_events[3]; cl_int err;
	err = clEnqueueReadBuffer( gpu_dev.commandQueue, bpd.d_rmsd_table,  CL_FALSE, 0, sizeof( float ) * bpd.n_good_node_pairs, &bpd.rms_data[1],   0, NULL, &read_buffer_events[0] );
	if ( err != CL_SUCCESS ) { utility_exit_with_message( "Failed to read from d_rmsd_table" ); }
	err = clEnqueueReadBuffer( gpu_dev.commandQueue, bpd.d_clash_table, CL_FALSE, 0, sizeof( float ) * bpd.n_good_node_pairs, &bpd.clash_data[1], 0, NULL, &read_buffer_events[1] );
	if ( err != CL_SUCCESS ) { utility_exit_with_message( "Failed to read from d_clash_table" ); }
	err = clEnqueueReadBuffer( gpu_dev.commandQueue, bpd.d_good_node_pairs, CL_FALSE, 0, sizeof( int ) * bpd.n_good_node_pairs, &bpd.good_node_pairs[1], 0, NULL, &read_buffer_events[2] );
	if ( err != CL_SUCCESS ) { utility_exit_with_message( "Failed to read from d_good_node_pairs" ); }

	clWaitForEvents( 3, read_buffer_events );

	// Don't delete the memory on the GPU: keep reusing it
	//gpu.Free( d_rmsd_table );
	//gpu.Free( d_clash_table );
	//gpu.Free( d_calculate_pair_table );

}

void
compare_gpu_result_against_cpu(
	helix_coord_map & helix_coords,
  bundle_helices_map & bundle_helices,
	helix_flipped_map & helix_flipped,
	utility::vector1< node_data > const & ndat,
	core::Size block1,
	core::Size block2,
  float maximum_rmsd_tolerance,
  float minimum_square_distance_tolerance,
	calculation_data_for_one_block_pair const & bpd
);


void
feed_rms_and_clash_scores_to_database(
	utility::sql_database::sessionOP db_session,
	utility::vector1< node_data > const & ndat,
	core::Size block1,
	core::Size block2,
	calculation_data_for_one_block_pair const & bpd // block-pair data
)
{
	std::string comparison_insert =
		"INSERT INTO node_comparisons(node_id_1, node_id_2, rmsd, clash_score) VALUES(?,?,?,?);";
	cppdb::statement comparison_insert_stmt(basic::database::safely_prepare_statement(comparison_insert,db_session));

  int block1_nnodes = (block1)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (block1-1)*BLOCK_SIZE;
  int block2_nnodes = (block2)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (block2-1)*BLOCK_SIZE;
  int block1_offset = (block1-1)*BLOCK_SIZE;
  int block2_offset = (block2-1)*BLOCK_SIZE;

	for ( int ii = 1; ii <= bpd.n_good_node_pairs; ++ii ) {
		int iipair = bpd.good_node_pairs[ ii ];
		if ( iipair / BLOCK_SIZE + 1 > block1_nnodes ) {
			std::cout << "Error iipair / BLOCK_SIZE + 1 > block1_nnodes: " << iipair << " " << iipair / BLOCK_SIZE << " " << block1_nnodes << std::endl;
			utility_exit_with_message( "bad node pair indices" );
		}
		if ( iipair % BLOCK_SIZE + 1 > block2_nnodes ) {
			std::cout << "Error iipair % BLOCK_SIZE + 1 > block2_nnodes: " << iipair << " " << iipair % BLOCK_SIZE << " " << block2_nnodes << std::endl;
			utility_exit_with_message( "bad node pair indices" );
		}
		int node1_index = iipair / BLOCK_SIZE + 1 + block1_offset;
		int node2_index = iipair % BLOCK_SIZE + 1 + block2_offset;
		//std::cout << "good pair " << ii << " n1: " << node1_index << " n2: " << node2_index << " (" << iipair /BLOCK_SIZE << "," << iipair % BLOCK_SIZE << ")" << std::endl;
		comparison_insert_stmt.bind(1, node1_index);
		comparison_insert_stmt.bind(2, node2_index);
		comparison_insert_stmt.bind(3, bpd.rms_data[ ii ] );
		comparison_insert_stmt.bind(4, bpd.clash_data[ ii ] );
		basic::database::safely_write_to_database(comparison_insert_stmt);
	}

	//db_session->begin();
	//for ( core::Size ii = 1; ii <= (core::Size) block1_nnodes; ++ii ) {
	//	core::Size iinode = ii + block1_offset;
	//	for ( core::Size jj = 1; jj <= (core::Size) block2_nnodes; ++jj ) {
	//		core::Size jjnode = jj + block2_offset;
	//		if ( ! calculate_pair( jj, ii ) ) continue;
	//		if ( rms(jj,ii) > 1.0 ) continue; // only save the good ones
	//		if ( clash(jj,ii) < 4.0 ) continue; // only save non-colliding helix pairs
	//		comparison_insert_stmt.bind(1, iinode);
	//		comparison_insert_stmt.bind(2, jjnode);
	//		comparison_insert_stmt.bind(3, rms(jj,ii) );
	//		comparison_insert_stmt.bind(4, clash(jj,ii) );
	//		basic::database::safely_write_to_database(comparison_insert_stmt);
	//	}
	//}
	//db_session->commit();

}


int
main( int argc, char * argv [] )
{
	using cppdb::statement;
	using cppdb::result;
	using namespace basic::database::schema_generator;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::wum;

	using namespace core;

	using ObjexxFCL::FArray2D;
	using ObjexxFCL::FArray1D;

	option.add( BundlePairRmsdCalculator::total_fractions, "The number of fractions to split the selection into");
	option.add( BundlePairRmsdCalculator::fraction, "The fraction of results to run RMSD comparisons on");
//	option.add( BundlePairRmsdCalculator::subfraction_size, "Size of each subfraction");

	// initialize core
	devel::init(argc, argv);

	// Initialize DB
	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	/**Get fraction of comparisons to be run**/
	int fraction = option[BundlePairRmsdCalculator::fraction].def(1);
	int total_fractions = option[BundlePairRmsdCalculator::total_fractions].def(1000);
//	int subfraction_size = option[BundlePairRmsdCalculator::subfraction_size].def(1000000);

	std::string get_count_string =
	"SELECT max(node_id) FROM helix_graph_nodes;";
	cppdb::statement get_count=basic::database::safely_prepare_statement(get_count_string, db_session);
	cppdb::result count_res=basic::database::safely_read_from_database(get_count);

	core::Size total_rows=0;
	while(count_res.next())
	{
		count_res >> total_rows;
	}

	core::Size num_rows_per_fraction = total_rows/total_fractions+1;

	core::Size overall_min_id = num_rows_per_fraction*fraction - num_rows_per_fraction + 1;
	core::Size overall_max_id = num_rows_per_fraction*fraction;
	TR << "overall_min_id: " << overall_min_id << std::endl;
	TR << "overall_max_id: " << overall_max_id << std::endl;

//	core::Size total_sub_fractions=num_rows_per_fraction/subfraction_size+1;
//	core::Size min_id=overall_min_id;
//	core::Size max_id=min_id+subfraction_size;
//
//	TR << "num_rows_per_fraction: " << num_rows_per_fraction << std::endl;
//	TR << "total_sub_fractions: " << total_sub_fractions << std::endl;

	// Load in the coordinates from the database
	helix_coord_map helix_coords;
	bundle_helices_map bundle_helices;
	helix_flipped_map helix_flipped;
	load_coordinates_from_database( db_session, helix_coords, bundle_helices, helix_flipped );
	TR << "Done populating helix coords and bundle helices" << std::endl;

	//*****SQL statements******//

	std::string count_nodes =
		"SELECT\n"
		" count(*)\n"
		"FROM protein_graph_nodes;\n";
	cppdb::statement count_nodes_stmt =
		basic::database::safely_prepare_statement(count_nodes, db_session);
	cppdb::result count_nodes_res =
		basic::database::safely_read_from_database( count_nodes_stmt );

	core::Size n_nodes;
	if ( ! count_nodes_res.next() ) {
		utility_exit_with_message( "Failed to retrieve the number of nodes from the database" );
	}
	count_nodes_res >> n_nodes;
	TR << "Done counting nodes from the database.  n_nodes = " << n_nodes  << std::endl;

	std::string select_all_nodes =
		"SELECT\n"
		"	node_id, substruct_id, comparison_element_1, comparison_element_2\n"
		"FROM protein_graph_nodes;\n";
	cppdb::statement select_all_nodes_stmt=
		basic::database::safely_prepare_statement(select_all_nodes, db_session);

	cppdb::result select_all_nodes_res=basic::database::safely_read_from_database(select_all_nodes_stmt);

	//*****END SQL statement*****//

	TR << "Done selecting nodes from the database"  << std::endl;

	utility::vector1< node_data > nodes( n_nodes );
	for ( core::Size ii = 1; ii <= n_nodes; ++ii ) {
		if ( ! select_all_nodes_res.next() ) {
			std::cerr << "Retrieved fewer nodes from the select_all_nodes query than the n_nodes query:" << ii << " vs " << n_nodes << std::endl;
			utility_exit_with_message( "problem with BundlePairRMSDCalculatorGPU application" );
		}
		select_all_nodes_res >> nodes[ii].node_index >> nodes[ii].bundle_index >> nodes[ii].helix_1_ind >> nodes[ii].helix_2_ind;
	}

	// OK: we'll divvy up the work into blocks of 256 x 256
	// For each block we'll first figure out whether we've already sent the coordinates
	// for the set of nodes on the "x axis" and for the set of nodes on the "y axis",
	// ship them to the GPU if necessary, and launch a kernel on the GPU to translate
	// them to the origin.
	// Then we'll compute an array of booleans (chars) to represent which of the
	// node/node comparisons should be performed and which shouldn't.
	// Then we'll allocate space on the GPU to hold the result of the computation
	// and ship the array of booleans to the GPU.
	// Then we'll launch a second kernel to compute the RMS and clash scores
	// making sure that the first kernel has first finished.

	basic::gpu::GPU gpu;
	std::vector< std::string > programs = gpu_rmsd_programs();
	gpu.RegisterProgram( programs );

	int const n_node_blocks = ( n_nodes + BLOCK_SIZE - 1 ) / BLOCK_SIZE;
	utility::vector1< int > coords_already_shipped( n_node_blocks, 0 );
	utility::vector1< gpu_node_block_ptrs > block_ptrs( n_node_blocks );

	int const n_block_pairs = n_node_blocks * ( n_node_blocks - 1 ) / 2 + n_node_blocks;
	int count_block_pairs_completed = 0;

	//ObjexxFCL::FArray2D< float > rms1( BLOCK_SIZE, BLOCK_SIZE );
	//ObjexxFCL::FArray2D< float > clash1( BLOCK_SIZE, BLOCK_SIZE );
	//ObjexxFCL::FArray2D< unsigned char > calculate_pair1( BLOCK_SIZE, BLOCK_SIZE );
	//cl_mem d_rms1, d_clash1, d_calculate_pair1;
	//
	//ObjexxFCL::FArray2D< float > rms2( BLOCK_SIZE, BLOCK_SIZE );
	//ObjexxFCL::FArray2D< float > clash2( BLOCK_SIZE, BLOCK_SIZE );
	//ObjexxFCL::FArray2D< unsigned char > calculate_pair2( BLOCK_SIZE, BLOCK_SIZE );
	//cl_mem d_rms2, d_clash2, d_calculate_pair2;
	//
	//cl_event compute_rms_event;

	utility::vector1< calculation_data_for_one_block_pair > block_pairs_data(2);
	for ( int ii = 1; ii <= 2; ++ii ) {
		// debug
		block_pairs_data[ ii ].rms_data_table.dimension( BLOCK_SIZE, BLOCK_SIZE );
		block_pairs_data[ ii ].clash_data_table.dimension( BLOCK_SIZE, BLOCK_SIZE );
		// end debug
		block_pairs_data[ ii ].rms_data.resize( BLOCK_SIZE * BLOCK_SIZE );
		block_pairs_data[ ii ].clash_data.resize( BLOCK_SIZE * BLOCK_SIZE );
		block_pairs_data[ ii ].good_node_pairs.resize( BLOCK_SIZE * BLOCK_SIZE );
		block_pairs_data[ ii ].n_good_node_pairs = 0;
		block_pairs_data[ ii ].d_rmsd_table = 0;
		block_pairs_data[ ii ].d_clash_table = 0;
		block_pairs_data[ ii ].d_calculate_pair_table = 0;
		block_pairs_data[ ii ].d_good_node_pairs = 0;
		block_pairs_data[ ii ].d_n_good_node_pairs = 0;
	}

	db_session->begin(); // write to the database once, only.  At the end, commit a single transaction.

	int which_index = 1; // either 1 or 2 -- keep the GPU loaded by switching between working on block_pairs_data[1] and block_pairs_data[2]
	int other_index = 2; // either 2 or 1 -- the opposite of which_index
	int lastii, lastjj;

	float maximum_tolerated_rms = 1.0f; // should probably be turned into a flag
	float minimum_acceptible_square_distance = 4.0; // should probably be turned into a flag -- 2.0A distance between Calpha and any other atom in the other helix

	for ( int ii = 1; ii <= n_node_blocks; ++ii ) {
		for ( int jj = ii; jj <= n_node_blocks; ++jj ) {

			// 1. send coordinates if necessary
			if ( ! coords_already_shipped[ ii ] ) {
				ship_coordinates_to_gpu( gpu, helix_coords,  bundle_helices, nodes, ii, block_ptrs[ ii ] );
				coords_already_shipped[ ii ] = 1;
			}
			if ( ! coords_already_shipped[ jj ] ) {
				ship_coordinates_to_gpu( gpu, helix_coords,  bundle_helices, nodes, jj, block_ptrs[ jj ] );
				coords_already_shipped[ jj ] = 1;
			}

			// launch a new kernel
			compute_rms_and_clash_score_for_block_pair(
				gpu, nodes, ii, jj, block_ptrs[ ii ], block_ptrs[ jj ],
				maximum_tolerated_rms, minimum_acceptible_square_distance,
				block_pairs_data[ which_index ] );
			//std::cout << "submitting rms/clash job for " << ii << " " << jj << " with index " << which_index << " and event: " << block_pairs_data[ which_index ].read_n_good_node_pairs_event << std::endl;

			// swap
			int tmp = which_index;
			which_index = other_index;
			other_index = tmp;

			// if this is not the first round, then pull the data from the last round from the GPU and send it to the database
			if ( ii != 1 || jj != 1 ) {

				//std::cout << "waiting for results from " << lastii << " " << lastjj << " with index " << which_index << " event: " << block_pairs_data[ which_index ].read_n_good_node_pairs_event << std::endl;;

				// wait until the last compute_rms_and_clash_score kenel finishes
				clWaitForEvents( 1, &block_pairs_data[ which_index ].read_n_good_node_pairs_event );
				//std::cout << "n_good_pairs: " << block_pairs_data[ which_index ].n_good_node_pairs << std::endl;

				// ok, now read back the good helix-pair data from the GPU and append it to the database session
				if ( block_pairs_data[ which_index ].n_good_node_pairs != 0 ) {
					retrieve_rms_and_clash_scores_from_gpu( gpu, block_pairs_data[ which_index ] );
					feed_rms_and_clash_scores_to_database( db_session, nodes, lastii, lastjj, block_pairs_data[ which_index ] );
				}
				//compare_gpu_result_against_cpu( helix_coords, bundle_helices, helix_flipped, nodes, lastii, lastjj,
				//	maximum_tolerated_rms, minimum_acceptible_square_distance, block_pairs_data[ which_index ] );
			}
			lastii = ii; lastjj = jj;

			++count_block_pairs_completed;
			if ( count_block_pairs_completed % 100 == 0 ) {
				std::cout << "Completed " << count_block_pairs_completed << " of " << n_block_pairs;
				std::cout << " (" << (float) count_block_pairs_completed / n_block_pairs << ")" << std::endl;
			}
		}
	}

	// very last iteration.  "other_index" holds the position of the last work unit submitted.
	// Wait until the final compute_rms_and_clash_score kenel finishes before exiting.
	clWaitForEvents( 1, &block_pairs_data[ other_index ].read_n_good_node_pairs_event );

	if ( block_pairs_data[ other_index ].n_good_node_pairs != 0 ) {
		retrieve_rms_and_clash_scores_from_gpu( gpu, block_pairs_data[ other_index ] );
		feed_rms_and_clash_scores_to_database( db_session, nodes, lastii, lastjj, block_pairs_data[ other_index ] );
		//compare_gpu_result_against_cpu( helix_coords, bundle_helices, helix_flipped, nodes, lastii, lastjj,
		//	maximum_tolerated_rms, minimum_acceptible_square_distance, block_pairs_data[ other_index ] );
	}


	db_session->commit(); // almost done: now commit all of the stored RMSD/clash-score pairs to the database

	return 0;
}

void
compare_gpu_result_against_cpu(
	helix_coord_map & helix_coords,
  bundle_helices_map & bundle_helices,
	helix_flipped_map & helix_flipped,
	utility::vector1< node_data > const & ndat,
	core::Size block1,
	core::Size block2,
  float maximum_rmsd_tolerance,
  float minimum_square_distance_tolerance,
	calculation_data_for_one_block_pair const & bpd
)
{
	using namespace ObjexxFCL;
	using core::Size;

  int block1_nnodes = (block1)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (block1-1)*BLOCK_SIZE;
  int block2_nnodes = (block2)*BLOCK_SIZE <= ndat.size() ? BLOCK_SIZE : ndat.size() - (block2-1)*BLOCK_SIZE;
  int block1_offset = (block1-1)*BLOCK_SIZE;
  int block2_offset = (block2-1)*BLOCK_SIZE;

	int n_wrong_clash(0), n_wrong_rms(0), n_correct_clash(0), n_correct_rms(0);

	std::cout << "compare gpu result against cpu: " << block1_nnodes << " " << block2_nnodes << std::endl;

	int count_cpu_good_pairs = 0;
	int count_gpu_good_pairs = 0;
	for ( int ii = 1; ii <= block1_nnodes; ++ii ) {
		for ( int jj = 1; jj <= block2_nnodes; ++jj ) {
			core::Size node_1_node_id, node_1_bundle_id, node_1_helix_id_1, node_1_helix_id_2;
			core::Size node_2_node_id, node_2_bundle_id, node_2_helix_id_1, node_2_helix_id_2;

			int iijj_indpair = (ii-1)*BLOCK_SIZE + (jj-1);

			//comparison_res >>
			//node_1_node_id >> node_1_bundle_id >> node_1_helix_id_1 >> node_1_helix_id_2 >>
			//node_2_node_id >> node_2_bundle_id >> node_2_helix_id_1 >> node_2_helix_id_2;

			// hijack the code Tim wrote:
			node_1_node_id = ii+block1_offset; node_2_node_id = jj+block2_offset;
			node_1_bundle_id = ndat[ node_1_node_id ].bundle_index; node_2_bundle_id = ndat[ node_2_node_id ].bundle_index;
			node_1_helix_id_1 = ndat[ node_1_node_id ].helix_1_ind; node_2_helix_id_1 = ndat[ node_2_node_id ].helix_1_ind;
			node_1_helix_id_2 = ndat[ node_1_node_id ].helix_2_ind; node_2_helix_id_2 = ndat[ node_2_node_id ].helix_2_ind;

			// skip if node_1_node_id >= node_2_node_id or if node1 and node2 are from the same bundle.
			//if ( node_1_node_id >= node_2_node_id ) continue;
			if ( node_1_bundle_id >= node_2_bundle_id ) {
				if ( count_gpu_good_pairs < bpd.n_good_node_pairs && bpd.good_node_pairs[ count_gpu_good_pairs+1 ] == iijj_indpair ) {
					std::cout << "oops! skipping over gpu good pair " << ii << " " << jj << " " << node_1_bundle_id << " " << node_2_bundle_id << "\n";
					std::cout << bpd.good_node_pairs[ count_gpu_good_pairs+1 ] << " " << bpd.good_node_pairs[ count_gpu_good_pairs+1 ] / BLOCK_SIZE;
					std::cout << " " << bpd.good_node_pairs[ count_gpu_good_pairs+1 ] % BLOCK_SIZE << " " << count_gpu_good_pairs << " " << bpd.good_node_pairs[ count_gpu_good_pairs+1 ] << " " << iijj_indpair << std::endl;
				}

				continue;
			}

			utility::vector1< numeric::xyzVector<core::Real> > node_1_coords;
			node_1_coords.insert( node_1_coords.end(), helix_coords[node_1_helix_id_1].begin(), helix_coords[node_1_helix_id_1].end() );
			node_1_coords.insert( node_1_coords.end(), helix_coords[node_1_helix_id_2].begin(), helix_coords[node_1_helix_id_2].end() );

			utility::vector1< numeric::xyzVector<core::Real> > node_2_coords;
			node_2_coords.insert( node_2_coords.end(), helix_coords[node_2_helix_id_1].begin(), helix_coords[node_2_helix_id_1].end() );
			node_2_coords.insert( node_2_coords.end(), helix_coords[node_2_helix_id_2].begin(), helix_coords[node_2_helix_id_2].end() );

			// compute center of mass for helices 1 and 2 for use later.
			numeric::xyzVector< numeric::Real > b1com(0.0), b2com(0.0);
			for ( Size i = 1; i <= node_1_coords.size(); ++i )
			{
				b1com += node_1_coords[i];
				b2com += node_2_coords[i];
			}
			b1com /= node_1_coords.size();
			b2com /= node_2_coords.size();

			//TESTING
//			TR << "Nodes " << node_1_node_id << " " << node_2_node_id << std::endl;
//			TR << "Helices 1: " << node_1_helix_id_1 << " " << node_1_helix_id_2 << std::endl;
//			TR << "Helices 2: " << node_2_helix_id_1 << " " << node_2_helix_id_2 << std::endl;
//			TR << "helix 1 coords size : " << node_1_coords.size() << std::endl;
//			TR << "helix 2 coords size : " << node_2_coords.size() << std::endl;
//			core::Real test_rmsd = numeric::model_quality::calc_rms(node_1_coords,node_2_coords);

			//Convert to FArrays
			runtime_assert(node_1_coords.size() == node_2_coords.size());

			//Save coords to FArrays
			FArray2D< numeric::Real > p1_coords( 3, node_1_coords.size() );
			FArray2D< numeric::Real > p2_coords( 3, node_2_coords.size() );
			for ( Size i = 1; i <= node_1_coords.size(); ++i )
			{
				for ( Size k = 1; k <= 3; ++k )// k = X, Y and Z
				{
					p1_coords(k,i) = node_1_coords[i](k);
					p2_coords(k,i) = node_2_coords[i](k);
				}
			}

			//uu is rotational matrix that transforms p2coords onto p1coords in a way that minimizes rms. Do this for ONLY the first two helices
			FArray1D< numeric::Real > ww( node_1_coords.size(), 1.0 );//weight matrix, all 1 for my purposes
			FArray2D< numeric::Real > uu( 3, 3, 0.0 );//transformation matrix
			numeric::Real ctx;
			numeric::model_quality::findUU( p1_coords, p2_coords, ww, node_1_coords.size(), uu, ctx );

			float rms_from_fast_rms;
			numeric::model_quality::calc_rms_fast( rms_from_fast_rms, p1_coords, p2_coords, ww, node_1_coords.size(), ctx );

			//Fill arrays with all coords not involved in the RMSD calculation. These coordinates
			//are used to calculate a clash score.
			utility::vector1< numeric::xyzVector<core::Real> > node_1_other_helix_coords;
			utility::vector1< core::Size> node_1_helices = bundle_helices[node_1_bundle_id];
			for(core::Size i=1; i<=node_1_helices.size(); ++i)
			{
				if( node_1_helices[i] != node_1_helix_id_1 &&
					node_1_helices[i] != node_1_helix_id_2)
				{
					if(helix_flipped[node_1_helices[i]])
					{
						node_1_other_helix_coords.insert(node_1_other_helix_coords.end(),
							helix_coords[node_1_helices[i]].begin(), helix_coords[node_1_helices[i]].end() );
					}
					else
					{
						utility::vector1< numeric::xyzVector<core::Real> > temp = helix_coords[node_1_helices[i]];
						reverse(temp.begin(), temp.end());
						node_1_other_helix_coords.insert( node_1_other_helix_coords.end(), temp.begin(), temp.end() );
					}
				}
			}

			utility::vector1< numeric::xyzVector<core::Real> > node_2_other_helix_coords;
			utility::vector1< core::Size> node_2_helices = bundle_helices[node_2_bundle_id];
			for(core::Size i=1; i<=node_2_helices.size(); ++i)
			{
				if( node_2_helices[i] != node_2_helix_id_1 &&
					node_2_helices[i] != node_2_helix_id_2)
				{
					if(helix_flipped[node_2_helices[i]])
					{
						node_2_other_helix_coords.insert(node_2_other_helix_coords.end(),
							helix_coords[node_2_helices[i]].begin(), helix_coords[node_2_helices[i]].end() );
					}
					else
					{
						utility::vector1< numeric::xyzVector<core::Real> > temp = helix_coords[node_2_helices[i]];
						reverse(temp.begin(), temp.end());
						node_2_other_helix_coords.insert( node_2_other_helix_coords.end(), temp.begin(), temp.end() );
					}
				}
			}

			assert(node_1_other_helix_coords.size() == node_2_other_helix_coords.size());

			//Add other helix coords to node vectors, these are now essentiall bundle-vectors
			node_1_coords.insert( node_1_coords.end(), node_1_other_helix_coords.begin(), node_1_other_helix_coords.end() );
			node_2_coords.insert( node_2_coords.end(), node_2_other_helix_coords.begin(), node_2_other_helix_coords.end() );

			//Save coords, now with third helix, to new FArrays
			FArray2D< numeric::Real > b1_coords( 3, node_1_coords.size() );
			FArray2D< numeric::Real > b2_coords( 3, node_2_coords.size() );
			for ( Size i = 1; i <= node_1_coords.size(); ++i )
			{
				for ( Size k = 1; k <= 3; ++k )// k = X, Y and Z
				{
					b1_coords(k,i) = node_1_coords[i](k);
					b2_coords(k,i) = node_2_coords[i](k);
				}
			}

			//move bundle to the origin using the center of mass for only the first two helices
			for ( int k = 1; k <= 3; ++k )
			{
				numeric::Real bundle_1_offset = 0.0;
				numeric::Real bundle_2_offset = 0.0;

				for ( int j = 1; j <= (int) (node_2_coords.size()-node_1_other_helix_coords.size()); ++j )
				{
					bundle_1_offset += b1_coords(k,j);
					bundle_2_offset += b2_coords(k,j);
				}
				bundle_1_offset /= (node_2_coords.size()-node_1_other_helix_coords.size());
				bundle_2_offset /= (node_2_coords.size()-node_1_other_helix_coords.size());

				for ( int j = 1; j <= (int) node_1_coords.size(); ++j )
				{
					b1_coords(k,j) -= bundle_1_offset;
					b2_coords(k,j) -= bundle_2_offset;
				}
			}

			//transform the coords of the entire bundle using rotational matrix produced by findUU for the first two helices
			FArray2D< numeric::Real > b2_coords_transformed(3, node_2_coords.size());
			for ( int i = 1; i <= (int) node_2_coords.size(); ++i )
			{
				b2_coords_transformed(1,i) = ( uu(1,1)*b2_coords(1,i) )+( uu(1,2)*b2_coords(2,i) ) +( uu(1,3)*b2_coords(3,i) );
				b2_coords_transformed(2,i) = ( uu(2,1)*b2_coords(1,i) )+( uu(2,2)*b2_coords(2,i) ) +( uu(2,3)*b2_coords(3,i) );
				b2_coords_transformed(3,i) = ( uu(3,1)*b2_coords(1,i) )+( uu(3,2)*b2_coords(2,i) ) +( uu(3,3)*b2_coords(3,i) );
			}

			//Calculate RMSD for the first two helices
			numeric::Real tot = 0;
			for ( int i = 1; i <= (int) (node_2_coords.size()-node_2_other_helix_coords.size()); ++i )
			{
				for ( int j = 1; j <= 3; ++j )
				{
					tot += std::pow( b1_coords(j,i) - b2_coords_transformed(j,i), 2 );
				}
			}
			core::Real cpu_bundle_pair_rmsd = std::sqrt(tot/(node_2_coords.size()-node_2_other_helix_coords.size()));

			////calc rms between transformed bundle_2_helix_3 pts and bundle_1_helix_3 pts
			//tot = 0;
			//for(int i = helix_coords[node_1_helix_id_1].size()+helix_coords[node_1_helix_id_2].size()+1;
			//	i <= node_2_coords.size(); ++i )
			//{
			//	for ( int j = 1; j <= 3; ++j ) {
			//		tot += std::pow( b1_coords(j,i) - b2_coords_transformed(j,i), 2 );
			//	}
			//}
			//core::Real third_helix_rmsd = std::sqrt(tot/node_1_other_helix_coords.size() );

			/// APL: calculate minimum distance between the calpha atoms of the middle residues of each of the
			/// helices and all atoms of the other helices for bundle 1 against 2, and for bundle 2 against 1.
			core::Real cpu_closest_contact = 12345;
			{ // scope
				utility::vector1< core::Size > const & iihelices = bundle_helices[ node_1_bundle_id ];
				utility::vector1< core::Size > const & jjhelices = bundle_helices[ node_2_bundle_id ];
				// first the center coordinates of the non-rms helices of node 1 against all atoms of the non-rms helices of node 2
				for ( core::Size kk = 1; kk <= iihelices.size(); ++kk ) {
					core::Size kkhelix_index = iihelices[kk];
					if ( kkhelix_index == node_1_helix_id_1 || kkhelix_index == node_1_helix_id_2 ) continue; // skip the RMS helices
					utility::vector1< numeric::xyzVector< core::Real> > const & kkhelix = helix_coords[ kkhelix_index ];
					core::Size ca_representative = (( kkhelix.size() / 4 ) / 2) * 4 + 1 + 1;
					numeric::xyzVector< core::Real > const kkcarep = kkhelix[ ca_representative ] - b1com;
					// now iterate across all the coordinates for the non-rms-helices in node 2
					for ( core::Size ll = 1; ll <= jjhelices.size(); ++ll ) {
						core::Size llhelix_index = jjhelices[ ll ];
						if ( llhelix_index == node_2_helix_id_1 || llhelix_index == node_2_helix_id_2 ) continue; // skip the RMS helices
						utility::vector1< numeric::xyzVector< core::Real> > const & llhelix = helix_coords[ llhelix_index ];
						for ( core::Size mm = 1; mm <= llhelix.size(); ++mm ) {
							numeric::xyzVector< core::Real > mm_coord = llhelix[mm] - b2com;
							numeric::xyzVector< core::Real > mm_transformed;
							mm_transformed(1) = ( uu(1,1)*mm_coord(1) )+( uu(1,2)*mm_coord(2) ) +( uu(1,3)*mm_coord(3) );
							mm_transformed(2) = ( uu(2,1)*mm_coord(1) )+( uu(2,2)*mm_coord(2) ) +( uu(2,3)*mm_coord(3) );
							mm_transformed(3) = ( uu(3,1)*mm_coord(1) )+( uu(3,2)*mm_coord(2) ) +( uu(3,3)*mm_coord(3) );
							core::Real d2 = kkcarep.distance_squared( mm_transformed  );
							cpu_closest_contact = d2 < cpu_closest_contact ? d2 : cpu_closest_contact;
						}
					}
				}
				// second the center coordinates of the non-rms helices of node 2 against all atoms of the non-rms helices of node 1
				for ( core::Size kk = 1; kk <= jjhelices.size(); ++kk ) {
					core::Size kkhelix_index = jjhelices[kk];
					if ( kkhelix_index == node_2_helix_id_1 || kkhelix_index == node_2_helix_id_2 ) continue; // skip the RMS helices
					utility::vector1< numeric::xyzVector< core::Real> > const & kkhelix = helix_coords[ kkhelix_index ];
					core::Size ca_representative = (( kkhelix.size() / 4 ) / 2) * 4 + 1 + 1;
					numeric::xyzVector< core::Real > kkcarep = kkhelix[ ca_representative ] - b2com;
					numeric::xyzVector< core::Real > kkcarep_transformed;
					kkcarep_transformed(1) = ( uu(1,1)*kkcarep(1) )+( uu(1,2)*kkcarep(2) ) +( uu(1,3)*kkcarep(3) );
					kkcarep_transformed(2) = ( uu(2,1)*kkcarep(1) )+( uu(2,2)*kkcarep(2) ) +( uu(2,3)*kkcarep(3) );
					kkcarep_transformed(3) = ( uu(3,1)*kkcarep(1) )+( uu(3,2)*kkcarep(2) ) +( uu(3,3)*kkcarep(3) );

					// now iterate across all the coordinates for the non-rms-helices in node 1
					for ( core::Size ll = 1; ll <= iihelices.size(); ++ll ) {
						core::Size llhelix_index = iihelices[ ll ];
						if ( llhelix_index == node_1_helix_id_1 || llhelix_index == node_1_helix_id_2 ) continue; // skip the RMS helices
						utility::vector1< numeric::xyzVector< core::Real> > const & llhelix = helix_coords[ llhelix_index ];
						for ( core::Size mm = 1; mm <= llhelix.size(); ++mm ) {
							core::Real d2 = kkcarep_transformed.distance_squared( llhelix[mm] - b1com );
							cpu_closest_contact = d2 < cpu_closest_contact ? d2 : cpu_closest_contact;
						}
					}
				}
			}

			//int iijj_indpair = (ii-1)*BLOCK_SIZE + (jj-1);
			bool iijjpair_good = cpu_bundle_pair_rmsd <= maximum_rmsd_tolerance && cpu_closest_contact >= minimum_square_distance_tolerance;

			// full debug -- look at all the values calculated on the GPU
			if ( std::abs( cpu_bundle_pair_rmsd - bpd.rms_data_table(jj,ii) ) > 1e-2 ) {
				std::cout << "pair: " << ii << " " << jj << " " << iijj_indpair << " rmsd calculated incorrectly on GPU: " << cpu_bundle_pair_rmsd << " " << bpd.rms_data_table(jj,ii) << " bundles: " << node_1_bundle_id << " " << node_2_bundle_id << std::endl;
			}
			if ( std::abs( cpu_closest_contact - bpd.clash_data_table(jj,ii) ) > 1e-2 ) {
				std::cout << "pair: " << ii << " " << jj << " " << iijj_indpair << " clash calculated incorrectly on GPU: " << cpu_closest_contact << " " << bpd.clash_data_table(jj,ii) << " bundles: " << node_1_bundle_id << " " << node_2_bundle_id << std::endl;
			}


			if ( iijjpair_good ) {
				++count_cpu_good_pairs;
				if ( count_gpu_good_pairs < bpd.n_good_node_pairs ) {
					if ( bpd.good_node_pairs[ count_gpu_good_pairs+1 ] == iijj_indpair ) {
						// ok -- we have a match between the CPU and the GPU
						++count_gpu_good_pairs;
						if ( std::abs( cpu_bundle_pair_rmsd - bpd.rms_data[ count_gpu_good_pairs ] ) > 1e-2 ) {
							++n_wrong_rms;
							std::cout << "wrong rmsd: " << iijj_indpair << " " << ii << " " << jj << " rmsd cpu: " << cpu_bundle_pair_rmsd << " count_gpu_good " << count_gpu_good_pairs << " rmsd gpu: " << bpd.rms_data[ count_gpu_good_pairs ] << std::endl;
						} else {
							++n_correct_rms;
						}

						if ( std::abs( cpu_closest_contact - bpd.clash_data[ count_gpu_good_pairs ] ) > 1e-3 ) {
							++n_wrong_clash;
							std::cout << "wrong clash: " << iijj_indpair << " " << ii << " " << jj << " clash cpu: " << cpu_closest_contact << " count_gpu_good " << count_gpu_good_pairs << " clash gpu: " << bpd.clash_data[ count_gpu_good_pairs ] << std::endl;
							//std::cout << "closest contact cpu: " << closest_contact << " gpu: " << clash_score_from_gpu( jj, ii ) << std::endl;
						} else {
							++n_correct_clash;
						}
					} else {
						// the GPU did not consider this to be a good pair, so hold it against the GPU so long as the CPU calculated
						// values are not near the threshold
            if ( std::abs( cpu_bundle_pair_rmsd - maximum_rmsd_tolerance ) > 1e-2 && std::abs( cpu_closest_contact - minimum_square_distance_tolerance ) > 1e-2 ) {
              ++n_wrong_rms; ++n_wrong_clash;
							std::cout << "gpu did not consider this pair good: " << iijj_indpair << " " << ii << " " << jj << " rms: " << cpu_bundle_pair_rmsd << " clash: " << cpu_closest_contact << std::endl;
            } else {
							// ok, give the GPU the benefit of the doubt even though we don't know what the GPU calculated
							++n_correct_rms; ++n_correct_clash;
						}
					}
				} else {
					// the GPU did not consier this to be a good pair
					// the GPU did not consider this to be a good pair, so hold it against the GPU so long as the CPU calculated
					// values are not near the threshold
					if ( std::abs( cpu_bundle_pair_rmsd - maximum_rmsd_tolerance ) > 1e-2 && std::abs( cpu_closest_contact - minimum_square_distance_tolerance ) > 1e-2 ) {
						++n_wrong_rms; ++n_wrong_clash;
						std::cout << "gpu did not consider this pair good: " << iijj_indpair << " " << ii << " " << jj << " rms: " << cpu_bundle_pair_rmsd << " clash: " << cpu_closest_contact << std::endl;
					} else {
						// ok, give the GPU the benefit of the doubt even though we don't know what the GPU calculated
						++n_correct_rms; ++n_correct_clash;
					}

				}
			} else {
				// bad pair: the GPU should also agree this is not a good pair
				if ( count_gpu_good_pairs < bpd.n_good_node_pairs ) {
					if ( bpd.good_node_pairs[ count_gpu_good_pairs+1 ] != iijj_indpair ) {
						// even though we don't know what the GPU calculated for this pair, give it the
						// benefit of the doubt
						++n_correct_rms; ++n_correct_clash;
					} else {
						// ok, it thinks the next pair is a good one, so, advance count_gpu_good_pairs, but
						// if the GPU and CPU calculated values disagree, then hold it against the GPU

            ++count_gpu_good_pairs;
            if ( std::abs( cpu_bundle_pair_rmsd - bpd.rms_data[ count_gpu_good_pairs ] ) > 1e-2 ) {
              ++n_wrong_rms;
							std::cout << "gpu wrongly considered this pair good: " << iijj_indpair << " " << ii << " " << jj << " rms: " << cpu_bundle_pair_rmsd << " clash: " << cpu_closest_contact << " and on the GPU: " << bpd.rms_data[ count_gpu_good_pairs ] << " " << bpd.clash_data[ count_gpu_good_pairs ] << std::endl;
            } else {
              ++n_correct_rms;
            }

            if ( std::abs( cpu_closest_contact - bpd.clash_data[ count_gpu_good_pairs ] ) > 1e-3 ) {
              ++n_wrong_clash;
							std::cout << "gpu wrongly considered this pair good: " << iijj_indpair << " " << ii << " " << jj << " rms: " << cpu_bundle_pair_rmsd << " clash: " << cpu_closest_contact << " and on the GPU: " << bpd.rms_data[ count_gpu_good_pairs ] << " " << bpd.clash_data[ count_gpu_good_pairs ] << std::endl;
            } else {
              ++n_correct_clash;
            }
					}
				} else {
					// there are no remaining good pairs on the GPU, so, give the GPU credit for not calling this a good
					// pair
					++n_correct_rms; ++n_correct_clash;
				}
			}
		}
	}

	if ( count_gpu_good_pairs != bpd.n_good_node_pairs ) {
		std::cout << "counting logic inside compare_gpu_result_against_cpu is incorrect: " << count_gpu_good_pairs << " != " << bpd.n_good_node_pairs << std::endl;
	}

	int count = n_wrong_clash + n_correct_clash;
	std::cout << "n comparisons: " << count << " rms: " << n_correct_rms << " vs " << n_wrong_rms << " ngoodcpu: " << count_cpu_good_pairs << std::endl;
	std::cout << "n comparisons: " << count << " clash: " << n_correct_clash << " vs " << n_wrong_clash << " ngoodcpu: " << count_cpu_good_pairs << std::endl;

}
