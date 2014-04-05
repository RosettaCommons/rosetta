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
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/chemical/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/stream_util.hh>
#include <basic/datacache/BasicDataCache.hh>
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
																Size const & res_to_add,
																Size const offset = 1 ){

	utility::vector1< Size > res_list_new(  res_list.size() + 1, 0 );

	for ( Size n = 1; n <= res_list.size(); n++ ){
		Size const m = res_list[ n ];
		if ( n < res_to_add )  res_list_new[ n ] = m;
		if ( n >= res_to_add ) res_list_new[ n+1 ] = m;
	}
	res_list_new[ res_to_add ] = res_list[ res_to_add ] - offset;

	return  res_list_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
reorder_res_list_after_append( utility::vector1< Size > const & res_list,
															 Size const & res_to_add,
															 Size const offset = 1 ){

	utility::vector1< Size > res_list_new(  res_list.size() + 1, 0 );

	for ( Size n = 1; n <= res_list.size(); n++ ){
		Size const m = res_list[ n ];
		if ( n < res_to_add )  res_list_new[ n ] = m;
		if ( n >= res_to_add ) res_list_new[ n+1 ] = m;
	}
	res_list_new[ res_to_add ] = res_list[ res_to_add-1 ] + offset;

	return  res_list_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_res_list_in_full_model_info_and_pdb_info( pose::Pose & pose, utility::vector1< Size > const & res_list_new ){
	FullModelInfoOP full_model_info = const_full_model_info( pose ).clone_info();
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
reorder_full_model_info_after_append( pose::Pose & pose, Size const res_to_add, Size const offset /* = 1 */ ){
	utility::vector1< Size > res_list_new = reorder_res_list_after_append( get_res_list_from_full_model_info_const( pose ), res_to_add, offset );
	update_res_list_in_full_model_info_and_pdb_info( pose, res_list_new );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
reorder_full_model_info_after_prepend( pose::Pose & pose, Size const res_to_add, Size const offset /* = 1 */ ){

	utility::vector1< Size > res_list_new = reorder_res_list_after_prepend( get_res_list_from_full_model_info_const( pose ), res_to_add, offset );
	update_res_list_in_full_model_info_and_pdb_info( pose, res_list_new );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_pdb_info_from_full_model_info( pose::Pose & pose ){

	using namespace core::pose;

	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );

	PDBInfoOP pdb_info = new PDBInfo( pose );
	pdb_info->set_numbering( res_list );

	// fill chains...
	utility::vector1< Size > chain_numbers = figure_out_chains_from_full_model_info( pose );
	utility::vector1< char > chains;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) chains.push_back( core::chemical::chr_chains[ chain_numbers[n]-1 ] );
	pdb_info->set_chains( chains );

	pose.pdb_info( pdb_info );

}

///////////////////////////////////////////////////////////////////////////////utility::vector1< Size >
utility::vector1< Size >
figure_out_chains_from_full_model_info( pose::Pose & pose ) {
	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
	return figure_out_chains_from_full_model_info_const( pose );
}


///////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_chains_full( pose::Pose const & pose ){

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	std::string const & sequence = full_model_info.full_sequence();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > cutpoint_open_in_full_model_including_terminus = cutpoint_open_in_full_model;
	if ( sequence.size() > 0 ) cutpoint_open_in_full_model_including_terminus.push_back( sequence.size() );
	runtime_assert( sequence.size() >= pose.total_residue() );

	Size start_res( 1 );
	utility::vector1< Size > chains_full;
	for ( Size i = 1; i <= cutpoint_open_in_full_model_including_terminus.size(); i++ ){
		Size const & end_res = cutpoint_open_in_full_model_including_terminus[ i ];
		for ( Size n = start_res; n <= end_res; n++ ) chains_full.push_back( i );
		start_res = end_res + 1;
	}
	return chains_full;

}

///////////////////////////////////////////////////////////////////////////////
// assign to different chains. [Note: we could actually use PDBInfo for this bookkeeping... perhaps that's a good idea]
// label them 1, 2, 3, etc.
utility::vector1< Size >
figure_out_chains_from_full_model_info_const( pose::Pose const & pose ) {

	using namespace core::pose::full_model_info;

	// first assign chains to full model
	utility::vector1< Size > const chains_full = get_chains_full( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );

	// now figure out chains in working model
	utility::vector1< Size > chains( pose.total_residue(), 1 );
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		if ( n > res_list.size() ) continue;
		runtime_assert( res_list[ n ] <= chains_full.size() );
		chains[ n ]  = chains_full[ res_list[ n ] ];
	}

	return chains;
}

///////////////////////////////////////////////////////////////////
bool
check_full_model_info_OK( pose::Pose const & pose ){

	using namespace core::pose::full_model_info;

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & conventional_numbering = full_model_info.conventional_numbering();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	std::string const & sequence = full_model_info.full_sequence();

	// very special case -- blank pose. could generalize to any pose with a virtual residue at end
	if ( res_list.size() == 0 && pose.total_residue() == 1 && pose.residue_type( 1 ).name3() == "XXX" ) return true;

	if ( res_list.size() != pose.total_residue() ) {
		TR << "res_list size != pose.total_residue() " << res_list.size() << " " << pose.total_residue() << std::endl;
		return false;
	}

	if ( sequence.size() != conventional_numbering.size() ) {
		TR << "sequence.size() != conventional_numbering.size()" << std::endl;
		return false;
	}


	if ( sequence.size() < pose.total_residue() ) {
		TR << "sequence.size() << pose.total_residue()" << std::endl;
		return false;
	}

	for ( Size n = 1; n <= res_list.size(); n++ ){

		Size const & res_num = res_list[ n ];
		runtime_assert( conventional_numbering.has_value( res_num ) );
		char sequence_char = sequence[ conventional_numbering.index( res_num )  - 1 ];

		if ( sequence_char == 'n' ) continue; // any nucleotide
		if ( sequence_char != pose.residue_type( n ).name1() ) {
			TR << "no match at " << n << " conventional numbering: " << res_num << "  sequence: " << sequence_char << " pose sequence: " <<  pose.residue_type( n ).name1() << std::endl;
			return false;
		}
	}

	return true;

}

	///////////////////////////////////////////////////////////////////////////////////////
	// mapping from working pose to full pose.
	utility::vector1< Size > const &
	get_res_list_from_full_model_info( pose::Pose & pose ) {
		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		return full_model_info.res_list();
	}

	///////////////////////////////////////////////////////////////////////////////////////
	// mapping from working pose to full pose.
	utility::vector1< Size > const &
	get_res_list_from_full_model_info_const( pose::Pose const & pose ) {
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		return full_model_info.res_list();
	}

	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< Size > >
	get_move_elements_from_full_model_info( pose::Pose & pose ){
		pose::full_model_info::make_sure_full_model_info_is_setup( pose );
		return get_move_elements_from_full_model_info_const( pose );
	}

	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< Size > >
	get_move_elements_from_full_model_info_const( pose::Pose const & pose ){

		FullModelInfo full_model_info = const_full_model_info( pose );
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

		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const & fixed_domain_map =  full_model_info.fixed_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();

		utility::vector1< Size > moving_res;
		for ( Size i = 1; i <= res_list.size(); i++ ) {
			if ( fixed_domain_map[ res_list[i] ] == 0 ) moving_res.push_back( i );
		}

		return moving_res;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_fixed_domain_from_full_model_info_const( pose::Pose const & pose ) {

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & fixed_domain_map =  full_model_info.fixed_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();

		utility::vector1< Size > fixed_domain_local;
		for ( Size i = 1; i <= res_list.size(); i++ ) {
			fixed_domain_local.push_back ( fixed_domain_map[ res_list[i] ] );
		}
		return fixed_domain_local;

	}

	/////////////////////////////////////////////////////////
	core::Size
	sub_to_full( core::Size const i, core::pose::Pose const & pose ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		runtime_assert( i > 0 && i <= full_model_info.res_list().size() );
		return full_model_info.res_list()[ i ];
	}

	/////////////////////////////////////////////////////////
	core::Size
	full_to_sub( core::Size const i, core::pose::Pose const & pose ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		return full_model_info.res_list().index( i );
	}

	/////////////////////////////////////////////////////////
	Size
	full_model_size( pose::Pose & pose ){
		return nonconst_full_model_info( pose ).size();
	}

	/////////////////////////////////////////////////////////
	void
	update_pose_domain_map( Pose & pose,
													Size & pose_domain_number,
													utility::vector1< Size > & pose_domain_map ){

		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		for ( Size k = 1; k <= res_list.size(); k++ ) pose_domain_map[ res_list[k] ] = pose_domain_number;

		utility::vector1< PoseOP > const & other_pose_list = full_model_info.other_pose_list();
		for ( Size n = 1; n <= other_pose_list.size(); n++ ){
			update_pose_domain_map( *(other_pose_list[ n ]), ++pose_domain_number, pose_domain_map );
		}

	}

	/////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_pose_domain_map( pose::Pose & pose ){

		utility::vector1< Size > pose_domain_map( full_model_size( pose ), 0 );

		Size pose_domain_number = 1;
		update_pose_domain_map( pose, pose_domain_number, pose_domain_map );
		return pose_domain_map;
	}


	/////////////////////////////////////////////////////////
	core::conformation::Residue const &
	get_residue( Size const seqpos_in_full_model,
							 pose::Pose const & pose,
							 bool & found_residue ){

		FullModelInfo const & full_model_info = const_full_model_info( pose );

		found_residue = false;

		utility::vector1< Size > const & res_list = full_model_info.res_list();
		if ( res_list.has_value( seqpos_in_full_model ) ){

			found_residue = true;
			return pose.residue( res_list.index( seqpos_in_full_model ) );

		} else {

			utility::vector1< PoseOP > const & other_pose_list = full_model_info.other_pose_list();
			for ( Size n = 1; n <= other_pose_list.size(); n++ ){
				core::conformation::Residue const & rsd = get_residue( seqpos_in_full_model, *(other_pose_list[ n ]), found_residue );
				if ( found_residue ) return rsd;
			}

		}

		return pose.residue( 1 ); // dummy return.
	}

	/////////////////////////////////////////////////////////
	core::conformation::Residue const &
	get_residue( Size const seqpos_in_full_model,
							 pose::Pose const & pose ){
		bool found_residue( false );
		core::conformation::Residue const & rsd = get_residue( seqpos_in_full_model, pose, found_residue );
		runtime_assert( found_residue );
		return rsd;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_res_num_from_pdb_info( pose::Pose const & pose ) {

		utility::vector1< Size > resnum;

		PDBInfoCOP pdb_info = pose.pdb_info();

		if ( pdb_info )	{
			for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( pdb_info->number( n ) );
		} else {
			for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( n );
		}

		return resnum;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	get_chain_for_full_model_resnum( Size const & resnum, pose::Pose const & pose ){
		return get_chains_full( pose )[ resnum ];
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	get_chain_for_resnum( Size const & resnum, pose::Pose const & pose ){
		return figure_out_chains_from_full_model_info_const( pose )[ resnum ];
	}

}
}
}
