// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/swa/FullModelInfo.cc
/// @brief  Mapping from a working pose into a bigger pose, for swa monte carlo stuff.
/// @author Rhiju Das

// Unit headers
#include <core/pose/full_model_info/FullModelInfoUtil.hh>

// Package headers

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

// C++
#include <string>
#include <map>

static basic::Tracer TR("core.pose.full_model_info.FullModelInfoUtil");

///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////

using namespace core;
using namespace core::pose::datacache;
using namespace basic::datacache;

namespace core {
namespace pose {
namespace full_model_info {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Following are utils and should not be in here...
utility::vector1<Size>
reorder_moving_res_list_after_delete( utility::vector1<Size> const & moving_res_list,
																			Size const & res_to_delete ){

	utility::vector1< Size > moving_res_list_new;

	for (Size i = 1; i <= moving_res_list.size(); i++ ){
		Size const n = moving_res_list[ i ];
		if ( n < res_to_delete ) moving_res_list_new.push_back( n );
		else if (n > res_to_delete ) moving_res_list_new.push_back( n-1 );
	}

	return moving_res_list_new;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
reorder_res_list_after_delete( utility::vector1< Size > const & res_list,
																	Size const & res_to_delete ){

	utility::vector1< Size > res_list_new(  res_list.size() - 1, 0 );

	for ( Size n = 1; n <= res_list.size(); n++ ){
		Size const m = res_list[ n ];
		if ( n < res_to_delete ) res_list_new[ n ] = m;
		else if ( n > res_to_delete ) res_list_new[ n-1 ] = m;
	}

	return res_list_new;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1<Size>
reorder_moving_res_list_after_insert( utility::vector1<Size> const & moving_res_list,
											 Size const & res_to_add ){

	utility::vector1< Size > moving_res_list_new;

	for (Size i = 1; i <= moving_res_list.size(); i++ ){
		Size const n = moving_res_list[ i ];
		if ( n < res_to_add ) moving_res_list_new.push_back( n );
	}
	moving_res_list_new.push_back( res_to_add );
	for (Size i = 1; i <= moving_res_list.size(); i++ ){
		Size const n = moving_res_list[ i ];
		if ( n >= res_to_add ) moving_res_list_new.push_back( n+1 );
	}

	return moving_res_list_new;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
reorder_res_list_after_prepend( utility::vector1< Size > const & res_list,
											Size const & res_to_add ){

	utility::vector1< Size > res_list_new(  res_list.size() + 1, 0 );

	for ( Size n = 1; n <= res_list.size(); n++ ){
		Size const m = res_list[ n ];
		if ( n < res_to_add )  res_list_new[ n ] = m;
		if ( n >= res_to_add ) res_list_new[ n+1 ] = m;
	}
	res_list_new[ res_to_add ] = res_list[ res_to_add ] - 1;

	return  res_list_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
reorder_res_list_after_append( utility::vector1< Size > const & res_list,
											Size const & res_to_add ){

	utility::vector1< Size > res_list_new(  res_list.size() + 1, 0 );

	for ( Size n = 1; n <= res_list.size(); n++ ){
		Size const m = res_list[ n ];
		if ( n < res_to_add )  res_list_new[ n ] = m;
		if ( n >= res_to_add ) res_list_new[ n+1 ] = m;
	}
	res_list_new[ res_to_add ] = res_list[ res_to_add-1 ]+1;

	return  res_list_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_res_list_in_full_model_info_and_pdb_info( pose::Pose & pose, utility::vector1< Size > const & res_list_new ){
	FullModelInfoOP full_model_info = const_full_model_info_from_pose( pose ).clone_info();
	full_model_info->set_res_list( res_list_new );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
	runtime_assert( check_full_model_info_OK( pose ) );
	update_pdb_info_from_full_model_info( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_delete( pose::Pose & pose, Size const res_to_delete ){

	utility::vector1< Size > res_list_new = reorder_res_list_after_delete( get_res_list_from_full_model_info_const( pose ), res_to_delete );
	update_res_list_in_full_model_info_and_pdb_info( pose, res_list_new );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_append( pose::Pose & pose, Size const res_to_add ){

	utility::vector1< Size > res_list_new = reorder_res_list_after_append( get_res_list_from_full_model_info_const( pose ), res_to_add );
	update_res_list_in_full_model_info_and_pdb_info( pose, res_list_new );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_prepend( pose::Pose & pose, Size const res_to_add ){

	utility::vector1< Size > res_list_new = reorder_res_list_after_prepend( get_res_list_from_full_model_info_const( pose ), res_to_add );
	update_res_list_in_full_model_info_and_pdb_info( pose, res_list_new );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_pdb_info_from_full_model_info( pose::Pose & pose ){

	using namespace core::pose;

	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );

	//	for ( Size i = 1; i <= res_list.size(); i++ ) std::cout << "SUB_TO_FULL " << res_list[i] << std::endl;

	//	PDBInfoOP pdb_info = pose.pdb_info();
	//	if ( ! pdb_info ) pdb_info = new PDBInfo( pose );
	PDBInfoOP pdb_info = new PDBInfo( pose );
	pdb_info->set_numbering( res_list );

	// fill chains...
	utility::vector1< Size > chain_numbers = figure_out_chains_from_full_model_info( pose );
	utility::vector1< char > chains;
	std::string const chain_char = "ABCDEFGHIJKLMNOPQRSTUVWZYX"; // hope this works.
	for ( Size n = 1; n <= pose.total_residue(); n++ ) chains.push_back( chain_char[ chain_numbers[n]-1 ] );
	pdb_info->set_chains( chains );

	pose.pdb_info( pdb_info );

}

///////////////////////////////////////////////////////////////////////////////
// assign to different chains. [Note: we could actually use PDBInfo for this bookkeeping... perhaps that's a good idea]
// label them 1, 2, 3, etc.
utility::vector1< Size >
figure_out_chains_from_full_model_info( pose::Pose & pose ) {

	using namespace core::pose::full_model_info;

	// first assign chains to full model
	utility::vector1< Size > chains_full;

	FullModelInfo const & full_model_info = nonconst_full_model_info_from_pose( pose );
	std::string const & sequence = full_model_info.full_sequence();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );

	utility::vector1< Size > cutpoint_open_in_full_model_including_terminus = cutpoint_open_in_full_model;
	cutpoint_open_in_full_model_including_terminus.push_back( sequence.size() );
	runtime_assert( sequence.size() >= pose.total_residue() );

	Size start_res( 1 );

	for ( Size i = 1; i <= cutpoint_open_in_full_model_including_terminus.size(); i++ ){

		Size const & end_res = cutpoint_open_in_full_model_including_terminus[ i ];

		//		TR << "chain " << i << ": " << start_res << " to " << end_res << std::endl;
		for ( Size n = start_res; n <= end_res; n++ ) chains_full.push_back( i );

		start_res = end_res + 1;
	}


	// now figure out chains in working model
	utility::vector1< Size > chains( pose.total_residue(), 1 );

	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		runtime_assert( res_list[ n ] <= chains_full.size() );
		chains[ n ]  = chains_full[ res_list[ n ] ];
		TR.Debug << "Setting chain at " << n << " to " << chains[ n  ] << std::endl;
	}

	return chains;
}


///////////////////////////////////////////////////////////////////
bool
check_full_model_info_OK( pose::Pose const & pose ){

	using namespace core::pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	std::string const & sequence = full_model_info.full_sequence();

	if ( res_list.size() != pose.total_residue() ) {
		TR.Debug << "res_list size != pose.total_residue() " << res_list.size() << " " << pose.total_residue() << std::endl;
		return false;
	}

	if ( sequence.size() < pose.total_residue() ) {
		TR.Debug << "sequence.size() << pose.total_residue()" << std::endl;
		return false;
	}

	for ( Size n = 1; n <= res_list.size(); n++ ){
		char sequence_char = sequence[ res_list[ n ] - 1 ];
		if ( sequence_char == 'n' ) continue; // any nucleotide
		if ( sequence_char != pose.residue_type( n ).name1() ) {
			TR.Debug << "no match at " << n << " " << sequence_char << ' ' <<  pose.residue_type( n ).name1() << std::endl;
			return false;
		}
	}

	return true;

}

	///////////////////////////////////////////////////////////////////////////////////////
	// mapping from working pose to full pose.
	utility::vector1< Size > const &
	get_res_list_from_full_model_info( pose::Pose & pose ) {
		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		return full_model_info.res_list();
	}

	///////////////////////////////////////////////////////////////////////////////////////
	// mapping from working pose to full pose.
	utility::vector1< Size > const &
	get_res_list_from_full_model_info_const( pose::Pose const & pose ) {
		FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
		return full_model_info.res_list();
	}


	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< Size > >
	get_move_elements_from_full_model_info( pose::Pose & pose ){


		FullModelInfo full_model_info = nonconst_full_model_info_from_pose( pose );
		utility::vector1< Size > const & fixed_domain_map =  full_model_info.fixed_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();

		utility::vector1< utility::vector1< Size > > move_elements;
		std::map< Size, utility::vector1< Size > > move_element_map;

		for ( Size i = 1; i <= res_list.size(); i++ ) {
			Size const i_full = res_list[ i ];
			Size const & domain = fixed_domain_map[ i_full ];
			if ( domain == 0 ) {
				// single residues
				move_elements.push_back( utility::tools::make_vector1( i_full ) );
			} else {
				// domains
				move_element_map[ domain ].push_back( i_full );
			}
		}

		for( std::map< Size, utility::vector1< Size > >::const_iterator iter = move_element_map.begin();
				 iter != move_element_map.end(); iter++ ){
			move_elements.push_back( iter->second );
		}

		return move_elements;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_moving_res_from_full_model_info( pose::Pose & pose ){

		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		utility::vector1< Size > const & fixed_domain_map =  full_model_info.fixed_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();

		utility::vector1< Size > moving_res;
		for ( Size i = 1; i <= res_list.size(); i++ ) {
			if ( fixed_domain_map[ res_list[i] ] == 0 ) moving_res.push_back( i );
		}

		return moving_res;
	}


	///////////////////////////////////////////////////////////////////////////////////////
	// I'm not happy with this -- copies code that is in PoseOP below...
	void
	fill_full_model_info_from_command_line( pose::Pose & pose ){

		using namespace core::sequence;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pose::full_model_info;

		if ( !option[ in::file::fasta ].user() ) return;

		// if pose is a subset of a bigger model, let's update that information.
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		std::string const target_sequence = fasta_sequence->sequence();

		utility::vector1< Size > input_res_list;
		if ( option[ in::file::input_res ].user() ) {

			input_res_list = option[ in::file::input_res ]();

			if ( input_res_list.size() != pose.total_residue() )  utility_exit_with_message( "Size of input_res does not match number of\
 residues in pose" );

			if ( input_res_list.size() > target_sequence.size() ) utility_exit_with_message( "Number in input_res exceeds number of resi\
dues in sequence" );

			for ( Size n = 1; n <= input_res_list.size(); n++ ){
				if ( target_sequence[ input_res_list[ n ] - 1 ] != pose.residue_type( n ).name1() ) {
					std::cerr << pose.sequence() << std::endl;
					std::cerr << target_sequence << std::endl;
					utility_exit_with_message( "pose sequence does not match fasta sequence, given input_res" );
				}
			}

		} else {
			for ( Size n = 1; n <= pose.total_residue(); n++ ) input_res_list.push_back( n );
		}

		utility::vector1< Size > cutpoint_open_in_full_model;
		if ( option[ full_model::cutpoint_open ].user() ) cutpoint_open_in_full_model = option[ full_model::cutpoint_open ]();

		FullModelInfoOP full_model_info =  new FullModelInfo( 	target_sequence );
		full_model_info->set_cutpoint_open_in_full_model( cutpoint_open_in_full_model );
		full_model_info->set_res_list( input_res_list );
		pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );

		// it might make sense to just use pdb_info instead of input_res_list...
		update_pdb_info_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.


	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	fill_full_model_info_from_command_line( utility::vector1< pose::PoseOP > input_poses ){

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( !option[ in::file::fasta ].user() ){
			for ( Size n = 1; n <= input_poses.size(); n++ ) nonconst_full_model_info_from_pose( *(input_poses[n]) );
			return;
		}

		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		std::string const desired_sequence = fasta_sequence->sequence();

		FullModelInfoOP full_model_info =	new FullModelInfo( desired_sequence );
		if( option[ full_model::cutpoint_open ].user() )	full_model_info->set_cutpoint_open_in_full_model( option[ full_model::cutpoint_open ]() );

		utility::vector1< Size > input_res_list = option[ in::file::input_res ]();
		utility::vector1< utility::vector1< Size > > pose_res_lists;
		utility::vector1< Size > domain_map( desired_sequence.size(), 0 );
		Size input_res_count( 0 );

		for ( Size n = 1; n <= input_poses.size(); n++ ) {

			Pose & pose = *(input_poses[n]);
			utility::vector1< Size > input_res_for_pose;

			for ( Size k = 1; k <= pose.total_residue(); k++ ){
				input_res_count++;
				runtime_assert( input_res_count <= input_res_list.size() );
				Size const & number_in_full_model = input_res_list[ input_res_count ];
				input_res_for_pose.push_back( number_in_full_model );
				domain_map[ number_in_full_model ] = n;
			}
			pose_res_lists.push_back( input_res_for_pose );
		}
		if ( input_res_count != input_res_list.size() ) utility_exit_with_message( "input_res size does not match pose size" );
		runtime_assert( input_res_count == input_res_list.size() );

		full_model_info->set_fixed_domain_map( domain_map );

		for ( Size n = 1; n <= input_poses.size(); n++ ) {
			Pose & pose = *(input_poses[n]);

			FullModelInfoOP full_model_info_for_pose = full_model_info->clone_info();
			full_model_info_for_pose->set_res_list( pose_res_lists[ n ] );

			// first pose holds information about other poses as its 'daughters' in a PoseTree.
			if ( n == 1 ){
				utility::vector1< PoseOP > other_pose_list;
				for ( Size n = 2; n <= input_poses.size(); n++ ) other_pose_list.push_back( input_poses[n] );
				full_model_info_for_pose->set_other_pose_list( other_pose_list );
			}

			pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info_for_pose );
			update_pdb_info_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.
		}


	}

	/////////////////////////////////////////////////////////
	core::Size
	sub_to_full( core::Size const i, core::pose::Pose & pose ){
		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		return full_model_info.res_list()[ i ];
	}

	/////////////////////////////////////////////////////////
	core::Size
	full_to_sub( core::Size const i, core::pose::Pose & pose ){
		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		return full_model_info.res_list().index( i );
	}


}
}
}
