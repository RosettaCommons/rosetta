// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Rosetta Headers
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh>
#include <protocols/farna/libraries/ChunkSet.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/base_pairs/BasePairStep.hh>
#include <protocols/farna/libraries/BasePairStepLibrary.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/toolbox/AtomID_Mapper.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/rna/util.hh>
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

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.libraries.RNA_ChunkLibrary" );

namespace protocols {
namespace farna {
namespace libraries {

using namespace core;
using namespace ObjexxFCL;
using namespace protocols::farna::base_pairs;

using core::Size;
using core::Real;

using core::pose::ResMap;

// magic number for special kind of chunk (base pair step) that is not user input,
//  but instead read in from Rosetta database.
Size const ROSETTA_LIBRARY_DOMAIN( 1000 );

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
RNA_ChunkLibrary::RNA_ChunkLibrary(){
	// currently nothing.
	chunk_coverage_ = 0.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
RNA_ChunkLibrary::RNA_ChunkLibrary( core::pose::Pose const & pose )
{
	utility::vector1< std::string > pdb_files_BLANK;
	utility::vector1< std::string > silent_files_BLANK;
	utility::vector1< core::Size > input_res_BLANK;
	initialize_rna_chunk_library( pdb_files_BLANK, silent_files_BLANK, pose, input_res_BLANK );
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
	utility::vector1< core::Size > const & input_res,
	utility::vector1< core::Size > const & allow_insert_res /* = blank */ )
{
	initialize_rna_chunk_library( pdb_files, silent_files, pose, input_res, allow_insert_res );
}

//////////////////////////////////////////////////////////////////////////////////////////////
// destructor
RNA_ChunkLibrary::~RNA_ChunkLibrary(){}

/// @brief clone the ChunkLibrary
RNA_ChunkLibraryOP
RNA_ChunkLibrary::clone() const
{
	RNA_ChunkLibraryOP chunk_library( new RNA_ChunkLibrary );
	*chunk_library = *this;
	return chunk_library;
}

//////////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_ChunkLibrary::num_chunks( Size const n ) const {
	return chunk_sets_[ n ]->num_chunks();
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::initialize_rna_chunk_library(
	utility::vector1 < std::string > const & pdb_files,
	utility::vector1 < std::string > const & silent_files,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & input_res,
	utility::vector1< core::Size > const & allow_insert_res /* = blank */ )
{
	std::string const & sequence_of_big_pose( pose.sequence() );
	coarse_rna_ = pose.residue( 1 ).is_coarse();
	do_rosetta_library_domain_check_ = true;

	// atom_level_domain_map keeps track of where chunks are placed -- only allow
	// fragment insertions *outside* these regions.
	if ( atom_level_domain_map_ == 0 ) {
		atom_level_domain_map_ = toolbox::AtomLevelDomainMapOP( new toolbox::AtomLevelDomainMap( pose, true /*map_to_vanilla*/, allow_insert_res ) );
	}
	covered_by_chunk_.dimension( sequence_of_big_pose.size(), false );

	utility::vector1< std::string > all_input_files;
	utility::vector1< bool > is_pdb_file;
	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		all_input_files.push_back( pdb_files[n] );
		is_pdb_file.push_back( true );
	}
	for ( Size n = 1; n <= silent_files.size(); n++ ) {
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
			if ( count > input_res.size() ) {
				std::cout << "Number of residues in scratch pose in RNA_ChunkLibrary " << count << " exceeds size of -input_res " << input_res.size() << std::endl;
				utility_exit_with_message( "problem with input_res" );
			}
			if ( input_res[ count ] < 1 ||
					input_res[ count ] > sequence_of_big_pose.size() ) {
				std::cout << "Problem with input_res: " << input_res[ count ] << " is bigger then length " << sequence_of_big_pose.size() << " of pose sequence " << sequence_of_big_pose << std::endl;
				utility_exit_with_message( "problem with input_res" );
			}
			if ( sequence_of_big_pose[ input_res[ count ] -1 ] != scratch_pose.sequence()[ i - 1 ] ) {
				std::cout << "Problem with input_file: " << all_input_files[n] << std::endl;
				std::cout << "mismatch in sequence   in  big pose: " << sequence_of_big_pose[ input_res[ count ] -1 ] << input_res[count] <<
					"  in input pose: " << scratch_pose.sequence()[ i - 1 ]  << i << std::endl;
				utility_exit_with_message( "mismatch in input_res sequence" );
			}
			res_map[ input_res[count ] ] = i;
		}

		ChunkSetOP chunk_set( new ChunkSet( pose_list, res_map ) );
		chunk_sets_.push_back( chunk_set );

		update_atom_level_domain_map( res_map, pose, scratch_pose, n );
		//check_fold_tree_OK( res_map, pose, scratch_pose );
	}

	if ( count != input_res.size() ) {
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
	chunk_sets_[ chunk_list_index ]->insert_chunk_into_pose( pose, chunk_pose_index, atom_level_domain_map_, do_rosetta_library_domain_check_ );
}


//////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
RNA_ChunkLibrary::get_indices_of_moving_chunks() const {
	utility::vector1< Size > moving_chunks;
	for ( Size n = 1; n <=  chunk_sets_.size(); n++ ) {
		if ( chunk_sets_[ n ]->num_chunks() > 1 )  moving_chunks.push_back( n );
	}
	return moving_chunks;
}

//////////////////////////////////////////////////////////////////////////////
Size
RNA_ChunkLibrary::num_moving_chunks() const { return get_indices_of_moving_chunks().size(); }


//////////////////////////////////////////////////////////////////////////////
bool
RNA_ChunkLibrary::random_chunk_insertion( core::pose::Pose & pose ) const
{

	utility::vector1< Size > const indices_of_moving_chunks = get_indices_of_moving_chunks();
	if ( indices_of_moving_chunks.size() == 0 ) return false;

	Size const chunk_set_index = numeric::random::rg().random_element( indices_of_moving_chunks );
	ChunkSet const & chunk_set( *chunk_sets_[ chunk_set_index ] );
	runtime_assert( chunk_set.num_chunks() > 1 );

	Size const chunk_index = static_cast <int> ( numeric::random::rg().uniform() * chunk_set.num_chunks() ) + 1;
	chunk_set.insert_chunk_into_pose( pose, chunk_index, atom_level_domain_map_, do_rosetta_library_domain_check_ );

	return true;
}

//////////////////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::update_atom_level_domain_map(
	ResMap const & res_map,
	core::pose::Pose const & pose,
	core::pose::Pose const & scratch_pose,
	core::Size const domain_num )
{
	using namespace core::id;
	using namespace core::conformation;

	// connected doesn't do anything anymore...
	FArray1D< bool > connected( pose.size(), false );
	covered_by_chunk_ = false;

	Size i_prev( 0 );

	for ( ResMap::const_iterator
			it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

		Size const i = it->first; //Index in big pose.
		Size const i_scratch = it->second; //Index in scratch pose (chunk).

		covered_by_chunk_( i ) = true;

		Residue const & rsd_i = pose.residue(i);
		for ( Size j = 1; j <= rsd_i.natoms(); j++ ) {

			std::string const & atomname = rsd_i.atom_name( j );
			Residue const & scratch_rsd = scratch_pose.residue(i_scratch);

			bool at_chainbreak =  ( i_scratch == 1 || scratch_pose.fold_tree().is_cutpoint( i_scratch - 1 ) || (i - i_prev ) > 1 );

			// awful, special case.
			if ( rsd_i.is_RNA() && involved_in_phosphate_torsion( atomname ) && at_chainbreak ) continue; // don't trust phosphates at beginning of chains.

			if ( scratch_rsd.has( atomname ) ) {
				Size const & scratch_index = scratch_pose.residue( i_scratch ).atom_index( atomname );
				AtomID const atom_id(j,i);
				if ( scratch_rsd.is_virtual( scratch_index ) ) continue;
				if ( !atom_level_domain_map_->has_domain( atom_id ) ) {
					// hack to get DNA/RNA compatibility -- has to do with using "vanilla" pose in atom_id_mapper as reference pose -- need to fix that.
					// why not just do name based lookup? -- rhiju, nov 2016.
				 	TR.Debug << TR.Cyan << "problem setting domain for " << rsd_i.name() << " " << rsd_i.atom_name( atom_id.atomno() ) << std::endl;
				 	continue;
				 }
				// special case: base pair steps should not overwrite user input domain
				if ( domain_num == ROSETTA_LIBRARY_DOMAIN &&
						atom_level_domain_map_->get_domain( atom_id ) > 0 ) continue;
				atom_level_domain_map_->set_domain( atom_id, domain_num);
			}
		}

		i_prev = i;
	}

}

//////////////////////////////////////////////////////////////////////////////
bool
RNA_ChunkLibrary::check_fold_tree_OK( pose::Pose const & pose ) const {

	for ( Size k = 1; k <= chunk_sets_.size(); k++ )  {
		ChunkSet & chunk_set = *(chunk_sets_[k]);
		bool const OK = chunk_set.check_fold_tree_OK( pose );
		if ( !OK ) {
			TR << TR.Red  << "Problem with pose fold tree -- not enough jumps to handle the number of chains in chunk set " << k << std::endl;
			TR << TR.Red  << "////////////////////////////////////////////////////////////////////" << std::endl;
			TR << TR.Red  << "WARNING!!! Ignoring this fold tree problem and continuing!!!" << std::endl;
			TR << TR.Red  << "////////////////////////////////////////////////////////////////////" << std::endl;
			//utility_exit_with_message( "FoldTree in pose does not have the right number of jumps to match chunk_res" );
		}
	}

	// Can't get here unless everything is OK!
	return true;
}

//////////////////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::figure_out_chunk_coverage()
{

	Size const tot_res( atom_level_domain_map_->atom_id_mapper()->nres() );
	runtime_assert( covered_by_chunk_.size() == tot_res );

	Size num_chunk_res( 0 );
	Size num_other_res( 0 );

	for ( Size n = 1; n <= tot_res; n++ ) {
		// Allow insert keeps track of where the chunk *aren't*, and
		// where other moves (fragments, jumps) can be carried out.
		if ( covered_by_chunk_(n) ) {
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
		if ( sequence[ i-1 ] != scratch_nt ) {
			std::cerr << "In full pose: " << i << " " << sequence[i-1] << "  In scratch pose: " << i_scratch_pose << " " << scratch_nt << std::endl;
			utility_exit_with_message(  "Mismatched sequence!!" );
			return false;
		}

	}
	return true;
}

////////////////////////////////////////////////////////////////
Size
RNA_ChunkLibrary::get_alignment_domain( pose::Pose const & pose ) const
{
	// figure out which domain might be good for alignment
	Size alignment_domain( 1 );
	Size anchor_rsd = get_anchor_rsd( pose ); // if there's a VRT residue, where it connects.
	if ( anchor_rsd > 0 ) {
		Size anchor_domain = atom_level_domain_map_->get_domain( anchor_rsd );
		if ( anchor_domain > 0 ) alignment_domain = anchor_domain;
	}
	//	TR << TR.Red << "alignment_domain: " << alignment_domain << " anchor_rsd " << anchor_rsd << std::endl;
	return alignment_domain;
}

////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::initialize_random_chunks( pose::Pose & pose, bool const dump_pdb /* = false */) const{

	if ( dump_pdb ) pose.dump_pdb( "start_"+string_of(0)+".pdb" );

	Size alignment_domain = get_alignment_domain( pose );

	for ( Size n = 1; n <= num_chunk_sets(); n++ ) {

		ChunkSet const & chunk_set( *chunk_sets_[ n ] );

		Size chunk_index = static_cast<int>( numeric::random::rg().uniform() * chunk_set.num_chunks() ) + 1;

		// JUST FOR TESTING
		if ( dump_pdb ) chunk_index = 1;

		//TR << "NUM_CHUNKS " << chunk_index << " " << chunk_set.num_chunks() << std::endl;
		chunk_set.insert_chunk_into_pose( pose, chunk_index, atom_level_domain_map_, do_rosetta_library_domain_check_ );

		// useful for tracking homology modeling: perhaps we can align to first chunk as well -- 3D alignment of Rosetta poses are
		// arbitrarily set to origin (except in special cases with virtual residues...)
		if ( n == alignment_domain && chunk_set.user_input()
				/*&&  pose.residue( pose.size() ).name3() != "VRT"*/ ) {
			align_to_chunk( pose, chunk_set, chunk_index  );
		}
		if ( dump_pdb ) pose.dump_pdb( "start_"+string_of(n)+".pdb" );
	}

}


///////////////////////////////////////////////////////////////////////
Size
RNA_ChunkLibrary::single_user_input_chunk() const
{
	Size user_input_chunk( 0 );
	for ( Size n = 1; n <= chunk_sets_.size(); n++ ) {
		ChunkSet const & chunk_set( *chunk_sets_[ n ] );
		if ( chunk_set.user_input() &&  chunk_set.num_chunks() == 1 ) {
			if ( user_input_chunk > 0 ) return 0; // more than one single user-input chunk
			user_input_chunk = n;
		}
	}
	return user_input_chunk;
}

///////////////////////////////////////////////////////////////////////
bool
RNA_ChunkLibrary::superimpose_to_single_user_input_chunk( pose::Pose & pose ) const{
	Size user_input_chunk = single_user_input_chunk();
	if ( user_input_chunk > 0 ) {
		align_to_chunk( pose, *(chunk_sets_[ user_input_chunk ]),  1  );
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::align_to_chunk( pose::Pose & pose, ChunkSet const & chunk_set, Size const chunk_index ) const{

	using namespace core::id;

	std::map< AtomID, AtomID > atom_id_map = chunk_set.get_atom_id_map( pose, *atom_level_domain_map_->atom_id_mapper() );

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
RNA_ChunkLibrary::set_atom_level_domain_map(toolbox::AtomLevelDomainMapOP atom_level_domain_map ){
	atom_level_domain_map_ = atom_level_domain_map;
}


//////////////////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::update_to_move_rosetta_library_chunks()
{
	atom_level_domain_map_->update_to_move_chunks_with_domain( libraries::ROSETTA_LIBRARY_DOMAIN );
	do_rosetta_library_domain_check_ = false; // now rosetta library positions are no longer marked with ROSETTA_LIBRARY_DOMAIN
}

//////////////////////////////////////////////////////////////////////////////
bool
check_base_pair_step_availability(
	BasePairStepLibrary const & base_pair_step_library,
	BasePairStepSequence const & base_pair_step_sequence )
{

	if ( base_pair_step_library.has_value( base_pair_step_sequence ) ) return true;

	std::string tag( base_pair_step_sequence.tag() );
	TR << TR.Red << base_pair_step_library.database_dir() << " does not have "  << base_pair_step_sequence.tag() << " for residue numbers " << tag;

	if ( !base_pair_step_library.canonical() ||
			std::count( tag.begin(), tag.end(), 'n' ) > 0 ) {
		TR << " so no base pair steps will be sampled there, just fragments." << std::endl;
		return false;
	}

	utility_exit_with_message( "Problem with canonical base pair step library." );
	return false;
}

//////////////////////////////////////////////////////////////////////////////
void
RNA_ChunkLibrary::setup_base_pair_step_chunks(
	pose::Pose const & pose,
	utility::vector1< BasePairStep > const & base_pair_steps,
	BasePairStepLibrary const & base_pair_step_library )
{
	using namespace core::id;
	using namespace core::pose;

	if ( chunk_sets_.size() > ( ROSETTA_LIBRARY_DOMAIN - 1 ) ) utility_exit_with_message( "hey update AtomLevelDomainMap & RNA_ChunkLibrary to not use magic numbers like 999 and 1000" );

	for ( Size m = 1; m <= base_pair_steps.size(); m++ ) {

		BasePairStep const & base_pair_step = base_pair_steps[ m ];

		BasePairStepSequence base_pair_step_sequence( pose.sequence(), base_pair_step );
		if ( !check_base_pair_step_availability( base_pair_step_library, base_pair_step_sequence ) ) continue;

		// happens if poses in the silent file in the database get filtered out.
		if ( base_pair_step_library.mini_pose_list( base_pair_step_sequence ).size() == 0 ) return;

		if ( !base_pair_step_moving( base_pair_step, atom_level_domain_map_, pose ) ) continue;

		ResMap res_map;
		res_map[ base_pair_step.i()      ] = 1;
		res_map[ base_pair_step.i_next() ] = 2;
		res_map[ base_pair_step.j()      ] = 3;
		res_map[ base_pair_step.j_next() ] = 4;

		ChunkSetOP chunk_set( new ChunkSet( base_pair_step_library.mini_pose_list( base_pair_step_sequence ),
			*base_pair_step_library.scratch_pose( base_pair_step_sequence ),
			res_map ) );
		chunk_set->set_user_input( false );
		chunk_sets_.push_back( chunk_set );

		pose::Pose const & scratch_pose = *( base_pair_step_library.scratch_pose( base_pair_step_sequence ) );
		check_res_map( res_map, scratch_pose, pose.sequence() );

		update_atom_level_domain_map( res_map, pose, scratch_pose, ROSETTA_LIBRARY_DOMAIN );
	}

	figure_out_chunk_coverage();

}


} //libraries
} //farna
} //protocols

