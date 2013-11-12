 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
 // vi: set ts=2 noet:
 //  CVS information:
 //  $Revision: 1.1.2.1 $
 //  $Date: 2005/11/07 21:05:35 $
 //  $Author: rhiju $
 // (c) Copyright Rosetta Commons Member Institutions.
 // (c) This file is part of the Rosetta software suite and is made available under license.
 // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
 // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
 // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/rna/RNA_ChunkLibrary.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/BasePairStep.hh>
#include <protocols/rna/BasePairStepLibrary.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>

static numeric::random::RandomGenerator RG(2380934);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.rna.RNA_ChunkLibrary" ) ;

namespace protocols{
namespace rna{

	using namespace core;
	using namespace ObjexxFCL;

	using core::Size;
	using core::Real;

	using core::pose::ResMap;

	///////////////////////////////////////////////////////////////////////
	ChunkSet::ChunkSet( utility::vector1< core::pose::MiniPoseOP > const & mini_pose_list,
											ResMap const & res_map ) {

		mini_pose_list_ = mini_pose_list;

		res_map_ = res_map;

		// not much information in mini_pose --> assume that all atoms are OK for copying.
		core::pose::MiniPose const & mini_pose = *(mini_pose_list[ 1 ]);
		for ( Size i = 1; i <= mini_pose.total_residue(); i++ ){
			for ( Size j = 1; j <= mini_pose.coords()[i].size(); j++ ){
				atom_id_mask_[ core::id::AtomID( j, i ) ] = true;
			}
		}

	}

	///////////////////////////////////////////////////////////////////////
	ChunkSet::ChunkSet( utility::vector1< core::pose::PoseOP > const & pose_list,
											ResMap const & res_map ) {

		for ( Size n = 1; n <= pose_list.size(); n++ ) {
			mini_pose_list_.push_back( core::pose::MiniPoseOP( new core::pose::MiniPose( *(pose_list[n]) ) ) );
		}
		res_map_ = res_map;

		core::pose::Pose const & pose = *( pose_list[1] );
		for ( Size i = 1; i <= pose.total_residue(); i++ ){

			core::conformation::Residue rsd = pose.residue( i );
			for ( Size j = 1; j <= rsd.natoms(); j++ ){
				atom_id_mask_[ core::id::AtomID( j, i ) ] = !rsd.is_virtual( j );
			}

			// special case for magnesium, which has a couple virtual atoms that need to get moved around and to define stubs.
			// not elegant, but I want to get this working.
			if ( rsd.name3() == " MG" ){
				for ( Size j = 1; j <= rsd.natoms(); j++ ) 	atom_id_mask_[ core::id::AtomID( j, i ) ] = true;
			}
		}


	}

	///////////////////////////////////////////////////////////////////////
	ChunkSet::~ChunkSet()	{}


	///////////////////////////////////////////////////////////////////////
	void
	ChunkSet::insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index,toolbox::AllowInsertOP const & allow_insert ) const{

		using namespace core::pose;
		using namespace core::id;

		core::pose::MiniPose const & scratch_pose ( *(mini_pose_list_[ chunk_pose_index ]) );

		//		TR << "SCRATCH_POSE " << scratch_pose.sequence() << ' ' << scratch_pose.fold_tree() << std::endl;

		std::map< AtomID, AtomID > atom_id_map = get_atom_id_map( pose, allow_insert );

		copy_dofs( pose, scratch_pose, atom_id_map  );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::map< id::AtomID, id::AtomID >
	ChunkSet::get_atom_id_map(  core::pose::Pose & pose, toolbox::AllowInsertOP const & allow_insert ) const{

		std::map< id::AtomID, id::AtomID > atom_id_map;

		allow_insert->calculate_atom_id_map( pose, res_map_, mini_pose_list_[1]->fold_tree(), atom_id_map );

		// This should prevent copying dofs for virtual phosphates, if they are tagged as such in the input silent files.
		filter_atom_id_map_with_mask( atom_id_map );

		return atom_id_map;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	ChunkSet::filter_atom_id_map_with_mask( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ) const{

		using namespace core::id;

		std::map< AtomID, AtomID > atom_id_map_new;

		for ( std::map< AtomID, AtomID >::const_iterator
						it=atom_id_map.begin(), it_end = atom_id_map.end(); it != it_end; ++it ) {

			AtomID const & insert_atom_id = it->first;
			AtomID const & source_atom_id = it->second;

			std::map< AtomID, bool >::const_iterator it_mask = atom_id_mask_.find( source_atom_id );
			if ( it_mask == atom_id_mask_.end() ) utility_exit_with_message( "Some problem with atom_id_mask in defining atom_id_map " );

			if ( !it_mask->second ) continue; // this source_atom_id is not allowed by mask, probably came from a virtual phosphate.

			atom_id_map_new[ insert_atom_id ] = source_atom_id;
		}

		atom_id_map = atom_id_map_new;

	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	core::pose::MiniPoseOP const
	ChunkSet::mini_pose( Size const idx ) const {
		return mini_pose_list_[ idx ];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	RNA_ChunkLibrary::RNA_ChunkLibrary(){
		// currently nothing.
		chunk_coverage_ = 0.0;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	RNA_ChunkLibrary::RNA_ChunkLibrary(
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res )
	{
		utility::vector1< std::string > pdb_files_BLANK;
		initialize_rna_chunk_library( pdb_files_BLANK, silent_files, pose, input_res );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	// constructor -- needs a list of silent files. Each silent file
	//  has solutions for a particular piece of the desired pose.
	RNA_ChunkLibrary::RNA_ChunkLibrary(
								utility::vector1 < std::string > const & pdb_files,
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res )
	{
		initialize_rna_chunk_library( pdb_files, silent_files, pose, input_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	// destructor
	RNA_ChunkLibrary::~RNA_ChunkLibrary(){}


	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::initialize_rna_chunk_library(
								utility::vector1 < std::string > const & pdb_files,
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res )
	{
		std::string const & sequence_of_big_pose( pose.sequence() );
		coarse_rna_ = pose.residue( 1 ).is_coarse();

		// allow_insert keeps track of where chunks are placed -- only allow
		// fragment insertions *outside* these regions.
		allow_insert_ = new toolbox::AllowInsert( pose );
		covered_by_chunk_.dimension( sequence_of_big_pose.size(), false );

		utility::vector1< std::string > all_input_files;
		utility::vector1< bool > is_pdb_file;
		for ( Size n = 1; n <= pdb_files.size(); n++ ){
			all_input_files.push_back( pdb_files[n] );
			is_pdb_file.push_back( true );
		}
		for ( Size n = 1; n <= silent_files.size(); n++ ){
			all_input_files.push_back( silent_files[n] );
			is_pdb_file.push_back( false );
		}

		Size count( 0 );
		for ( Size n = 1; n <= all_input_files.size(); n++ ) {

			utility::vector1< pose::PoseOP > pose_list;
			process_input_file( all_input_files[n], pose_list, is_pdb_file[n], coarse_rna_ );

			core::pose::Pose const & scratch_pose( *(pose_list[1]) );

			// There may be more than one part of the pose to which this sequence maps.
			ResMap res_map;

			for ( Size i = 1; i <= scratch_pose.sequence().size(); i++ ) {
				count++;
				if ( sequence_of_big_pose[ input_res[ count ] -1 ] != scratch_pose.sequence()[ i - 1 ] ){
					std::cout << "Problem with input_file: " << all_input_files[n] << std::endl;
					std::cout << "mismatch in sequence   in  big pose: " << sequence_of_big_pose[ input_res[ count ] -1 ] << input_res[count] <<
						"  in input pose: " << scratch_pose.sequence()[ i - 1 ]  << i << std::endl;
					utility_exit_with_message( "mismatch in input_res sequence" );
				}
				res_map[ input_res[count ] ] = i;
			}

			ChunkSetOP chunk_set( new ChunkSet( pose_list, res_map ) );
			chunk_sets_.push_back( chunk_set );

			update_allow_insert( res_map, pose, scratch_pose, n );
			//check_fold_tree_OK( res_map, pose, scratch_pose );
		}

		if ( count != input_res.size() ){
			utility_exit_with_message( "Number of input res does not match total res in input silent files!" );
		}

		figure_out_chunk_coverage();

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::add_chunk_set(
												std::string const & silent_file,
												ResMap const & res_map,
												pose::Pose const & big_pose )
	{
		utility::vector1< pose::PoseOP > pose_list;

		process_input_file( silent_file, pose_list, false /*is_pdb*/, coarse_rna_ );
		check_res_map( res_map, *(pose_list[1]), big_pose.sequence() );

		ChunkSetOP chunk_set( new ChunkSet( pose_list, res_map ) );
		chunk_sets_.push_back( chunk_set );
	}

	////////////////////////////////////////////////////////////////////////////
	void RNA_ChunkLibrary::insert_chunk_into_pose(
					 pose::Pose & pose,
					 Size const & chunk_list_index,
					 Size const & chunk_pose_index ) const
	{
		chunk_sets_[ chunk_list_index ]->insert_chunk_into_pose( pose, chunk_pose_index, allow_insert_ );
	}


	//////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	RNA_ChunkLibrary::get_indices_of_moving_chunks() const {
		utility::vector1< Size > moving_chunks;
		for ( Size n = 1; n <=  chunk_sets_.size(); n++ ){
			if ( chunk_sets_[ n ]->num_chunks() > 1 )  moving_chunks.push_back( n );
		}
		return moving_chunks;
	}

	//////////////////////////////////////////////////////////////////////////////
	Size
	RNA_ChunkLibrary::num_moving_chunks() const { return get_indices_of_moving_chunks().size(); }


	//////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::random_chunk_insertion( core::pose::Pose & pose ) const{

		utility::vector1< Size > const indices_of_moving_chunks = get_indices_of_moving_chunks();
		if ( indices_of_moving_chunks.size() == 0 ) return false;

		Size const chunk_set_index = RG.random_element( indices_of_moving_chunks );
		ChunkSet const & chunk_set( *chunk_sets_[ chunk_set_index ] );
		runtime_assert( chunk_set.num_chunks() > 1 );

		Size const chunk_index = static_cast <int> ( RG.uniform() * chunk_set.num_chunks() ) + 1;
		chunk_set.insert_chunk_into_pose( pose, chunk_index, allow_insert_ );

		return true;
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::update_allow_insert( ResMap const & res_map,
																					 core::pose::Pose const & pose,
																					 core::pose::Pose const & scratch_pose,
																					 core::Size const domain_num )
	{
		using namespace core::id;
		using namespace core::conformation;

		// connected doesn't do anything anymore...
		FArray1D< bool > connected( pose.total_residue(), false );
		covered_by_chunk_ = false;

		Size i_prev( 0 );

		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

			Size const i = it->first; //Index in big pose.
			Size const i_scratch = it->second; //Index in scratch pose (chunk).

			covered_by_chunk_( i ) = true;

			Residue const & rsd_i = pose.residue(i);
			for ( Size j = 1; j <= rsd_i.natoms(); j++ ){

				std::string const & atomname = rsd_i.atom_name( j );
				Residue const & scratch_rsd = scratch_pose.residue(i_scratch);

				bool at_chainbreak =  ( i_scratch == 1 || scratch_pose.fold_tree().is_cutpoint( i_scratch - 1 ) ||	(i - i_prev ) > 1 );

				if ( rsd_i.is_RNA() && involved_in_phosphate_torsion( atomname ) && at_chainbreak ) continue; // don't trust phosphates at beginning of chains.

				if ( scratch_rsd.has( atomname ) ) {
					Size const & scratch_index = scratch_pose.residue( i_scratch ).atom_index( atomname );
					if ( !scratch_rsd.is_virtual( scratch_index ) ) {
						allow_insert_->set_domain( AtomID(j,i), domain_num);
					}
				}
			}

			i_prev = i;
		}

	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::check_fold_tree_OK( pose::Pose const & pose ){

		for (Size k = 1; k <= chunk_sets_.size(); k++ )  {
			ChunkSet & chunk_set = *(chunk_sets_[k]);
			bool const OK = chunk_set.check_fold_tree_OK( pose );
			if (!OK){
				std::cout << "Problem with pose fold tree -- not enough jumps to handle the number of chains in chunk set " << k << std::endl;
				utility_exit_with_message( "FoldTree in pose does not have the right number of jumps to match chunk_res" );
			}
		}

		// Can't get here unless everything is OK!
		return true;
	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	ChunkSet::check_fold_tree_OK( pose::Pose const & pose ){

		// Check where the chunk is mapped to in the big pose.
		// There should be at least the same number of jumps in the big pose
		//  as there are chains in the scratch_pose.
		utility::vector1< bool > is_chunk_res( pose.total_residue(), false );
		for ( ResMap::const_iterator
						it=res_map_.begin(), it_end = res_map_.end(); it != it_end; ++it ) {
			Size const i = it->first; //Index in big pose.
			is_chunk_res[ i ] = true;
		}

		Size const num_jumps_scratch = mini_pose_list_[1]->fold_tree().num_jump(); // number of chains - 1

		Size num_jumps_in_big_pose_in_scratch_region( 0 );
		for ( Size n = 1; n <= pose.num_jump(); n++ ) {
			if (! is_chunk_res[ pose.fold_tree().upstream_jump_residue( n ) ] ) continue;
			if (! is_chunk_res[ pose.fold_tree().downstream_jump_residue( n ) ] ) continue;
			num_jumps_in_big_pose_in_scratch_region++;
		}

		if ( num_jumps_scratch > num_jumps_in_big_pose_in_scratch_region ){
			std::cout << "Number of jumps in chunk pose               : " << num_jumps_scratch << std::endl;
			std::cout << "Number of jumps in full pose in chunk region: " << num_jumps_in_big_pose_in_scratch_region  << "  out of total jumps " << pose.num_jump() << std::endl;
			return false;
		}

		if ( num_jumps_scratch < num_jumps_in_big_pose_in_scratch_region ){
			//			std::cout << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
			//			std::cout << "Number of jumps in chunk pose               : " << num_jumps_scratch << std::endl;
			//			std::cout << "Does not match:" << std::endl;
			//			std::cout << "Number of jumps in full pose in chunk region: " << num_jumps_in_big_pose_in_scratch_region  << "  out of total jumps " << pose.num_jump() << std::endl;
			//			std::cout << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
			// Just a warning
			//return false;
		}

		return true;

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::figure_out_chunk_coverage()
	{

		Size const tot_res( allow_insert_->nres() );
		Size num_chunk_res( 0 );
		Size num_other_res( 0 );

		for ( Size n = 1; n <= tot_res; n++ ) {
			// Allow insert keeps track of where the chunk *aren't*, and
			// where other moves (fragments, jumps) can be carried out.
			if ( covered_by_chunk_(n) ){
				num_chunk_res++;
			} else {
				num_other_res++;
			}
		}
		chunk_coverage_ = Real( 3 * num_chunk_res ) / ( 3 * num_chunk_res +  tot_res );

	}


	//////////////////////////////////////////////////////////////////////////////
	// Check sequence here -- exit if no match.
	bool
	RNA_ChunkLibrary::check_res_map( ResMap const & res_map, pose::Pose const & scratch_pose, std::string const & sequence ) const{

		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

			// For now, just do bonded atoms...update later to do jumps too.
			Size const i = it->first; //Index in big pose.
			Size const i_scratch_pose = it->second; // Index in the little "chunk" or "scratch" pose

			char scratch_nt = scratch_pose.residue( i_scratch_pose ).name1();
			if ( sequence[ i-1 ] != scratch_nt ){
				std::cerr << "In full pose: " << i << " " << sequence[i-1] << "  In scratch pose: " << i_scratch_pose << " " << scratch_nt << std::endl;
				utility_exit_with_message(  "Mismatched sequence!!" );
				return false;
			}

		}
		return true;
	}

	////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::initialize_random_chunks( pose::Pose & pose, bool const dump_pdb /* = false */) const{
		for ( Size n = 1; n <= num_chunk_sets(); n++ ) {

			ChunkSet const & chunk_set( *chunk_sets_[ n ] );

			Size chunk_index = static_cast<int>( RG.uniform() * chunk_set.num_chunks() ) + 1;

			// JUST FOR TESTING
			if ( dump_pdb ) chunk_index = 1;

			//TR << "NUM_CHUNKS " << chunk_index << " " << chunk_set.num_chunks() << std::endl;
			chunk_set.insert_chunk_into_pose( pose, chunk_index, allow_insert_ );

			// useful for tracking homology modeling: perhaps we can align to first chunk as well -- 3D alignment of Rosetta poses are
			// arbitrarily set to origin (except in special cases with virtual residues...)
			if ( n==1  /*&&  pose.residue( pose.total_residue() ).name3() != "VRT"*/ ) align_to_chunk( pose, chunk_set, chunk_index  );

			if ( dump_pdb ) pose.dump_pdb( "start_"+string_of(n)+".pdb" );

		}

		//exit( 0 );

	}

	///////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::superimpose_to_first_chunk( pose::Pose & pose ) const{
		runtime_assert( chunk_sets_.size() > 0 );
		ChunkSet const & chunk_set( *chunk_sets_[ 1 ] );
		align_to_chunk( pose, chunk_set,  1  );

	}

	///////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::align_to_chunk( pose::Pose & pose, ChunkSet const & chunk_set, Size const chunk_index ) const{

		using namespace core::id;

		std::map< AtomID, AtomID > atom_id_map = chunk_set.get_atom_id_map( pose, allow_insert_ );

		id::AtomID_Map< id::AtomID >  alignment_atom_id_map; // weird alternative format needed for superimpose_pose
		core::pose::initialize_atomid_map( alignment_atom_id_map, pose, id::BOGUS_ATOM_ID );
		for ( std::map< AtomID, AtomID >::const_iterator
						it=atom_id_map.begin(), it_end = atom_id_map.end(); it != it_end; ++it ) {
			alignment_atom_id_map.set( it->first, it->second );
		}

		core::scoring::superimpose_pose( pose, *(chunk_set.mini_pose( chunk_index )), alignment_atom_id_map );

	}



	////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::set_allow_insert(toolbox::AllowInsertOP allow_insert ){
		allow_insert_ = allow_insert;
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::setup_base_pair_step_chunks( pose::Pose const & pose, utility::vector1< BasePairStep > base_pair_steps ){

		using namespace core::id;
		using namespace core::pose;

		if ( !base_pair_step_library_ ) base_pair_step_library_ = new BasePairStepLibrary;
		base_pair_step_library_->initialize();

		toolbox::AllowInsertOP allow_insert_original = allow_insert_->clone();

		if ( chunk_sets_.size() > 900 ) utility_exit_with_message( "hey update AllowInsert & RNA_ChunkLibrary to not use 999 as a magic number" );
		Size q( 1000 ); // heh heh, magic number.

		for ( Size m = 1; m <= base_pair_steps.size(); m++ ){

			BasePairStep const & base_pair_step = base_pair_steps[ m ];

			BasePairStepSequence base_pair_step_sequence( pose.sequence(), base_pair_step );
			runtime_assert( base_pair_step_library_->has_value( base_pair_step_sequence ) );

			if ( !allow_insert_original->get( named_atom_id_to_atom_id( NamedAtomID( " C1'", base_pair_step.i() ), pose ) ) ) continue;
			if ( !allow_insert_original->get( named_atom_id_to_atom_id( NamedAtomID( " C1'", base_pair_step.i_next() ), pose ) ) ) continue;
			if ( !allow_insert_original->get( named_atom_id_to_atom_id( NamedAtomID( " C1'", base_pair_step.j() ), pose ) ) ) continue;
			if ( !allow_insert_original->get( named_atom_id_to_atom_id( NamedAtomID( " C1'", base_pair_step.j_next() ), pose  ) ) ) continue;

			ResMap res_map;
			res_map[ base_pair_step.i()      ] = 1;
			res_map[ base_pair_step.i_next() ] = 2;
			res_map[ base_pair_step.j()      ] = 3;
			res_map[ base_pair_step.j_next() ] = 4;

			ChunkSetOP chunk_set( new ChunkSet( base_pair_step_library_->mini_pose_list( base_pair_step_sequence ), res_map ) );
			chunk_sets_.push_back( chunk_set );

			pose::Pose const & scratch_pose = *base_pair_step_library_->scratch_pose( base_pair_step_sequence );
			check_res_map( res_map, scratch_pose, pose.sequence() );

			q++;
			update_allow_insert( res_map, pose, scratch_pose, q );
		}

		figure_out_chunk_coverage();

	}


}
}

