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

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

// C++
#include <string>
#include <map>

static basic::Tracer tr("core.pose.full_model_info.FullModelInfoUtil");

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
reorder_sub_to_full_after_delete( utility::vector1< Size > const & sub_to_full,
																	Size const & res_to_delete ){

	utility::vector1< Size > sub_to_full_new(  sub_to_full.size() - 1, 0 );

	for ( Size n = 1; n <= sub_to_full.size(); n++ ){
		Size const m = sub_to_full[ n ];
		if ( n < res_to_delete ) sub_to_full_new[ n ] = m;
		else if ( n > res_to_delete ) sub_to_full_new[ n-1 ] = m;
	}

	return sub_to_full_new;

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
reorder_sub_to_full_after_prepend( utility::vector1< Size > const & sub_to_full,
											Size const & res_to_add ){

	utility::vector1< Size > sub_to_full_new(  sub_to_full.size() + 1, 0 );

	for ( Size n = 1; n <= sub_to_full.size(); n++ ){
		Size const m = sub_to_full[ n ];
		if ( n < res_to_add )  sub_to_full_new[ n ] = m;
		if ( n >= res_to_add ) sub_to_full_new[ n+1 ] = m;
	}
	sub_to_full_new[ res_to_add ] = sub_to_full[ res_to_add ] - 1;

	return  sub_to_full_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
reorder_sub_to_full_after_append( utility::vector1< Size > const & sub_to_full,
											Size const & res_to_add ){

	utility::vector1< Size > sub_to_full_new(  sub_to_full.size() + 1, 0 );

	for ( Size n = 1; n <= sub_to_full.size(); n++ ){
		Size const m = sub_to_full[ n ];
		if ( n < res_to_add )  sub_to_full_new[ n ] = m;
		if ( n >= res_to_add ) sub_to_full_new[ n+1 ] = m;
	}
	sub_to_full_new[ res_to_add ] = sub_to_full[ res_to_add-1 ]+1;

	return  sub_to_full_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_delete( pose::Pose & pose, Size const res_to_delete ){

	FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );

	utility::vector1< Size > sub_to_full_new = reorder_sub_to_full_after_delete( full_model_info.sub_to_full(), res_to_delete );
	utility::vector1< Size > moving_res_list_new = reorder_moving_res_list_after_delete( full_model_info.moving_res_list(), res_to_delete );

	full_model_info.set_sub_to_full( sub_to_full_new );
	full_model_info.set_moving_res_list( moving_res_list_new );

	update_pdb_info_from_sub_to_full( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_append( pose::Pose & pose, Size const res_to_add ){

	FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );

	utility::vector1< Size > sub_to_full_new = reorder_sub_to_full_after_append( full_model_info.sub_to_full(), res_to_add );
	utility::vector1< Size > moving_res_list_new = reorder_moving_res_list_after_insert( full_model_info.moving_res_list(), res_to_add );

	full_model_info.set_sub_to_full( sub_to_full_new );
	full_model_info.set_moving_res_list( moving_res_list_new );

	update_pdb_info_from_sub_to_full( pose );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_prepend( pose::Pose & pose, Size const res_to_add ){

	FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );

	utility::vector1< Size > sub_to_full_new = reorder_sub_to_full_after_prepend( full_model_info.sub_to_full(), res_to_add );
	utility::vector1< Size > moving_res_list_new = reorder_moving_res_list_after_insert( full_model_info.moving_res_list(), res_to_add );

	full_model_info.set_sub_to_full( sub_to_full_new );
	full_model_info.set_moving_res_list( moving_res_list_new );

	update_pdb_info_from_sub_to_full( pose );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_pdb_info_from_sub_to_full( pose::Pose & pose ){

	using namespace core::pose;

	FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
	utility::vector1< Size > sub_to_full = full_model_info.sub_to_full();

	// following may not be necessary anymore -- sub_to_full *is* working_res
	utility::vector1< Size > working_res;
	for( Size n = 1; n <= pose.total_residue(); n++ ) working_res.push_back( sub_to_full[ n ] );

	PDBInfoOP pdb_info = pose.pdb_info();
	if ( ! pdb_info ) pdb_info = new PDBInfo( pose );

	pdb_info->set_numbering( working_res );
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
	utility::vector1< Size > const & sub_to_full = full_model_info.sub_to_full();
	std::string const & sequence = full_model_info.full_sequence();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();

	tr.Debug << "HEY! GOT FULL_MODEL_INFO " << sub_to_full.size() << " " << cutpoint_open_in_full_model.size() << std::endl;

	utility::vector1< Size > cutpoint_open_in_full_model_including_terminus = cutpoint_open_in_full_model;
	cutpoint_open_in_full_model_including_terminus.push_back( sequence.size() );
	runtime_assert( sequence.size() >= pose.total_residue() );

	Size start_res( 1 );

	for ( Size i = 1; i <= cutpoint_open_in_full_model_including_terminus.size(); i++ ){

		Size const & end_res = cutpoint_open_in_full_model_including_terminus[ i ];

		tr.Debug << "chain " << i << ": " << start_res << " to " << end_res << std::endl;
		for ( Size n = start_res; n <= end_res; n++ ) chains_full.push_back( i );

		start_res = end_res + 1;
	}


	// now figure out chains in working model
	utility::vector1< Size > chains( pose.total_residue(), 1 );

	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		chains[ n ]  = chains_full[ sub_to_full[ n ] ];
		tr.Debug << "Setting chain at " << n << " to " << chains[ n  ] << std::endl;
	}

	return chains;
}



//////////////////////////////////////////////////////////////////
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

		if ( input_res_list.size() != pose.total_residue() )	utility_exit_with_message( "Size of input_res does not match number of residues in pose" );

		if ( input_res_list.size() > target_sequence.size() ) utility_exit_with_message( "Number in input_res exceeds number of residues in sequence" );

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

	utility::vector1< Size > start_moving_res_list /*blank*/;

	utility::vector1< Size > cutpoint_open_in_full_model;
	if ( option[ full_model::cutpoint_open ].user() ) cutpoint_open_in_full_model = option[ full_model::cutpoint_open ]();

	FullModelInfoOP full_model_info_op =	new FullModelInfo(  input_res_list, start_moving_res_list, target_sequence, cutpoint_open_in_full_model );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info_op );

	// it might make sense to just use pdb_info instead of input_res_list...
	update_pdb_info_from_sub_to_full( pose ); // for output pdb or silent file -- residue numbering.


}





}
}
}
