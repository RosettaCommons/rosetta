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
#include <core/pose/full_model_info/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/types.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/loop_graph/LoopGraph.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/stream_util.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

// C++
#include <string>
#include <map>

static thread_local basic::Tracer TR( "core.pose.full_model_info.util" );

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
update_pdb_info_from_full_model_info( pose::Pose & pose ){

	using namespace core::pose;

	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );

	PDBInfoOP pdb_info( new PDBInfo( pose ) );
	pdb_info->set_numbering( const_full_model_info( pose ).full_model_parameters()->full_to_conventional( res_list ) );
	pdb_info->set_chains(    figure_out_conventional_chains_from_full_model_info( pose ) );

	pose.pdb_info( pdb_info );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_constraint_set_from_full_model_info( pose::Pose & pose ){

	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::id;

	FullModelParametersCOP full_model_parameters = const_full_model_info( pose ).full_model_parameters();
	ConstraintSetOP cst_set( new ConstraintSet );
	std::string const cst_string = full_model_parameters->cst_string();

	if ( cst_string.size() > 0 ) {
		full_model_parameters->update_pose_and_cst_set_from_cst_string( pose.residue( 1 ).residue_type_set() );
		Pose const & full_model_pose = full_model_parameters->full_model_pose_for_constraints();
		ConstraintSetCOP full_model_cst_set = full_model_parameters->cst_set();

		utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
		id::SequenceMappingOP sequence_map( new SequenceMapping );
		for ( Size n = 1; n <= full_model_pose.total_residue(); n++ ){
			if ( res_list.has_value( n ) ){
				sequence_map->push_back( res_list.index( n ) );
			} else {
				sequence_map->push_back( 0 );
			}
		}
		cst_set = full_model_cst_set->remapped_clone( full_model_pose, pose, sequence_map );
	}

	pose.constraint_set( cst_set );
}

/////////////////////////////////////////////////////////////////////
void
update_disulfides_from_full_model_info( pose::Pose & pose ){
	FullModelInfo const & full_model_info = nonconst_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	FullModelParametersCOP full_model_parameters = full_model_info.full_model_parameters();
	utility::vector1< std::pair< Size, Size > > const & disulf_bonds =
		full_model_parameters->get_res_list_as_pairs( DISULFIDE );
	utility::vector1< std::pair< Size, Size > > working_disulf_bonds;

	for ( Size i = 1; i <= disulf_bonds.size(); i++ ) {
		Size const & res1_full = disulf_bonds[ i ].first;
		Size const & res2_full = disulf_bonds[ i ].second;
		if ( !res_list.has_value( res1_full ) ) continue;
		if ( !res_list.has_value( res2_full ) ) continue;
		std::pair< Size, Size > const disulfide_pair =
			std::make_pair( full_model_info.full_to_sub( res1_full ),
											full_model_info.full_to_sub( res2_full ) );
		working_disulf_bonds.push_back( disulfide_pair );
		runtime_assert( pose.residue_type( disulfide_pair.first  ).is_sidechain_thiol() || pose.residue_type( disulfide_pair.first  ).is_disulfide_bonded() );
		runtime_assert( pose.residue_type( disulfide_pair.second ).is_sidechain_thiol() || pose.residue_type( disulfide_pair.second ).is_disulfide_bonded() );
	}

	pose.conformation().fix_disulfides( working_disulf_bonds );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_pose_objects_from_full_model_info( pose::Pose & pose ) {
	update_pdb_info_from_full_model_info( pose );
	update_constraint_set_from_full_model_info( pose );
	update_disulfides_from_full_model_info( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
update_res_list_in_full_model_info_and_pdb_info( pose::Pose & pose, utility::vector1< Size > const & res_list_new ){
	FullModelInfoOP full_model_info = const_full_model_info( pose ).clone_info();
	full_model_info->set_res_list( res_list_new );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
	runtime_assert( check_full_model_info_OK( pose ) );
	update_pose_objects_from_full_model_info( pose );
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


///////////////////////////////////////////////////////////////////////////////
utility::vector1< char >
figure_out_conventional_chains_from_full_model_info( pose::Pose const & pose ) {
	utility::vector1< char > chains;
	utility::vector1< char > const & conventional_chains = const_full_model_info( pose ).conventional_chains();
	if ( conventional_chains.size() > 0 ){
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		for ( Size n = 1; n <= pose.total_residue(); n++ ) chains.push_back( conventional_chains[ res_list[n] ] );
	} else {
		utility::vector1< Size > chain_numbers = figure_out_chain_numbers_from_full_model_info_const( pose );
		for ( Size n = 1; n <= pose.total_residue(); n++ ) chains.push_back( core::chemical::chr_chains[ chain_numbers[n]-1 ] );
	}
	return chains;
}

///////////////////////////////////////////////////////////////////////////////utility::vector1< Size >
utility::vector1< Size >
figure_out_chain_numbers_from_full_model_info( pose::Pose & pose ) {
	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
	return figure_out_chain_numbers_from_full_model_info_const( pose );
}


///////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_chains_full( pose::Pose const & pose ){

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	std::string const & sequence = full_model_info.full_sequence();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	return get_chains_from_cutpoint_open( cutpoint_open_in_full_model, sequence.size() );
}

///////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_chains_from_cutpoint_open( utility::vector1< Size > const & cutpoint_open, Size const nres ) {

	utility::vector1< Size > cutpoint_open_including_terminus = cutpoint_open;
	if ( nres > 0 ) cutpoint_open_including_terminus.push_back( nres );

	Size start_res( 1 );
	utility::vector1< Size > chains;
	for ( Size i = 1; i <= cutpoint_open_including_terminus.size(); i++ ){
		Size const & end_res = cutpoint_open_including_terminus[ i ];
		for ( Size n = start_res; n <= end_res; n++ ) chains.push_back( i );
		start_res = end_res + 1;
	}
	return chains;
}

///////////////////////////////////////////////////////////////////////////////
// assign to different chains. [Note: we could actually use PDBInfo for this bookkeeping... perhaps that's a good idea]
// label them 1, 2, 3, etc.
utility::vector1< Size >
figure_out_chain_numbers_from_full_model_info_const( pose::Pose const & pose ) {

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
utility::vector1< Size >
figure_out_dock_domain_map_from_full_model_info_const( pose::Pose const & pose ) {

	using namespace core::pose::full_model_info;

	// first assign chains to full model
	utility::vector1< Size > const dock_domain_map = const_full_model_info( pose ).dock_domain_map();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );

	// now figure out chains in working model
	utility::vector1< Size > dock_domain_map_for_pose;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) dock_domain_map_for_pose.push_back( dock_domain_map[ res_list[ n ] ] );

	return dock_domain_map_for_pose;
}

///////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_sample_res_for_pose( pose::Pose const & pose ) {
	FullModelInfo const & full_model_info( const_full_model_info( pose ) );
	utility::vector1< Size > const & sample_res = full_model_info.sample_res();
	return full_model_info.full_to_sub( sample_res );
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
		char sequence_char = sequence[ res_num   - 1 ];

		if ( sequence_char == 'n' ) continue; // any nucleotide
		if ( sequence_char != pose.residue_type( n ).name1() ) {
			TR << res_list << std::endl;
			TR << "no match at " << n << " conventional numbering: " << res_num << "  sequence: " << sequence_char << " pose sequence: " <<  pose.residue_type( n ).name1() << std::endl;
			TR << "POSE SEQUENCE " << pose.sequence() << std::endl;
			TR << "FULL MODEL SEQUENCE " << sequence << std::endl;
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
		runtime_assert( full_model_info_defined( pose ) );
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		return full_model_info.res_list();
	}

	// this function works even if full_model_info is not defined.
	utility::vector1< Size >
	get_res_list_const( pose::Pose const & pose ) {
		if ( full_model_info_defined( pose ) ) return get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > res_list;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) res_list.push_back( n );
		return res_list;
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
		utility::vector1< Size > const & input_domain_map =  full_model_info.input_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();

		utility::vector1< utility::vector1< Size > > move_elements;
		std::map< Size, utility::vector1< Size > > move_element_map;

		for ( Size i = 1; i <= res_list.size(); i++ ) {
			Size const i_full = res_list[ i ];
			Size const & domain = input_domain_map[ i_full ];
			if ( domain == 0 ) {
				// single residues
				move_elements.push_back( utility::tools::make_vector1( i_full ) );
			} else {
				// domains
				move_element_map[ domain ].push_back( i_full );
			}
		}

		for( std::map< Size, utility::vector1< Size > >::const_iterator iter = move_element_map.begin(), end = move_element_map.end();
				 iter != end; ++iter ){
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
	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_input_domain_from_full_model_info_const( pose::Pose const & pose ) {

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & input_domain_map =  full_model_info.input_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();

		utility::vector1< Size > input_domain_local;
		for ( Size i = 1; i <= res_list.size(); i++ ) {
			input_domain_local.push_back ( input_domain_map[ res_list[i] ] );
		}
		return input_domain_local;

	}

	/////////////////////////////////////////////////////////
	core::Size
	sub_to_full( core::Size const i, core::pose::Pose const & pose ){
		utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
		runtime_assert( i >= 1 && i <= res_list.size() );
		return res_list[ i ];
	}

	utility::vector1< Size >
	sub_to_full( utility::vector1< Size > const & res_list, core::pose::Pose const & pose ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		return full_model_info.sub_to_full( res_list );
	}

	/////////////////////////////////////////////////////////
	core::Size
	full_to_sub( core::Size const i, core::pose::Pose const & pose ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		return full_model_info.res_list().index( i );
	}

	/////////////////////////////////////////////////////////
	utility::vector1< Size >
	full_to_sub( utility::vector1< Size > const & res_list, core::pose::Pose const & pose ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		return full_model_info.full_to_sub( res_list );
	}

	/////////////////////////////////////////////////////////
	Size
	full_model_size( pose::Pose & pose ){
		return nonconst_full_model_info( pose ).size();
	}

	/////////////////////////////////////////////////////////
	void
	update_pose_domain_map_const( Pose const & pose,
																Size & pose_domain_number,
																utility::vector1< Size > & pose_domain_map ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		for ( Size k = 1; k <= res_list.size(); k++ ) {
			pose_domain_map[ res_list[k] ] = pose_domain_number;
		}
		utility::vector1< PoseOP > const & other_pose_list = full_model_info.other_pose_list();
		for ( Size n = 1; n <= other_pose_list.size(); n++ ){
			update_pose_domain_map_const( *(other_pose_list[ n ]), ++pose_domain_number, pose_domain_map );
		}
	}

	/////////////////////////////////////////////////////////
	void
	update_pose_domain_map( Pose & pose,
													Size & pose_domain_number,
													utility::vector1< Size > & pose_domain_map ){
		make_sure_full_model_info_is_setup( pose );
		update_pose_domain_map_const( pose, pose_domain_number, pose_domain_map );
	}


	/////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_pose_domain_map( pose::Pose & pose ){
		make_sure_full_model_info_is_setup( pose );
		return figure_out_pose_domain_map_const( pose );
	}

	/////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_pose_domain_map_const( pose::Pose const & pose ){

		utility::vector1< Size > pose_domain_map( const_full_model_info( pose ).size(), 0 );
		Size pose_domain_number = 1;
		update_pose_domain_map_const( pose, pose_domain_number, pose_domain_map );
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
	utility::vector1< int >
	get_res_num_from_pdb_info( pose::Pose const & pose ) {

		utility::vector1< int > resnum;

		PDBInfoCOP pdb_info = pose.pdb_info();

		if ( ( pdb_info != 0 ) && !pdb_info->obsolete())	{
			for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( pdb_info->number( n ) );
		} else {
			for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( n );
		}
		return resnum;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< char >
	get_chains_from_pdb_info( pose::Pose const & pose ) {

		utility::vector1< char > chains;

		PDBInfoCOP pdb_info = pose.pdb_info();

		if ( ( pdb_info != 0 ) && !pdb_info->obsolete())	{
			for ( Size n = 1; n <= pose.total_residue(); n++ ) chains.push_back( pdb_info->chain( n ) );
		} else {
			for ( Size n = 1; n <= pose.total_residue(); n++ ) chains.push_back( ' ' );
		}
		return chains;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	get_chain_for_full_model_resnum( Size const & resnum, pose::Pose const & pose ){
		return get_chains_full( pose )[ resnum ];
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	get_chain_for_resnum( Size const & resnum, pose::Pose const & pose ){
		return figure_out_chain_numbers_from_full_model_info_const( pose )[ resnum ];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// deprecate soon in favor of direct loop-graph call
	//
	/// @brief Finds the number of missing residues in a pose.
	/// @details This function returns the number of missing residues in the pose.
	/// The pose is passed by nonconst reference, so that the full_model_info can be
	/// setup, if needed.
	Size
	get_number_missing_residues_and_connections( pose::Pose & pose ) {
		make_sure_full_model_info_is_setup( pose );
		utility::vector1< char > missing_residues;
		utility::vector1< utility::vector1< Size > > loop_suite_lengths;
		return get_number_missing_residues_and_connections( pose, missing_residues, loop_suite_lengths );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// deprecate soon in favor of direct loop-graph call
	//
	/// @brief Finds the number of missing residues in a pose.
	/// @details This function returns the number of missing residues in the pose.
	/// The missing_residues vector is passed by nonconst reference, so its values
	/// can be modified and accessed by this function and the calling method.
	Size
	get_number_missing_residues_and_connections( pose::Pose const & pose,
																							 utility::vector1< char > & missing_residues ){
		utility::vector1< utility::vector1< Size > > loop_suites;
		return get_number_missing_residues_and_connections( pose, missing_residues, loop_suites );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// deprecate soon in favor of direct loop-graph call
	//
	/// @brief Finds the number of missing residues in a pose.
	/// @details This function returns the number of missing residues in the pose.
	/// The loop_suites vector is passed by nonconst reference, so its values
	/// can be modified and accessed by this function and the calling method.
	Size
	get_number_missing_residues_and_connections( pose::Pose const & pose,
																							 utility::vector1< utility::vector1< Size > > & loop_suites ){
		utility::vector1< char > missing_residues;
		return get_number_missing_residues_and_connections( pose, missing_residues, loop_suites );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// sort of a hodge-podge function -- might be better to unify with LoopGraph?
	//
	/// @brief Finds the number of missing residues in a pose.
	/// @details This function returns the number of missing residues in the pose.
	/// The missing_res and loop_suites vectors are passed by nonconst reference, so their
	/// values can be modified and accessed by this function and the calling method.
	Size
	get_number_missing_residues_and_connections( pose::Pose const & pose,
																							 utility::vector1< char > & missing_residues,
																							 utility::vector1< utility::vector1< Size > > & loop_suites ) {

		////////////////////////////////////////////////////////
		// following is old.
		// probably should replace with LoopGraph -- see below.
		// REMOVE in 2015.
		////////////////////////////////////////////////////////
		using namespace core::pose::full_model_info;
		utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map_const( pose );
		utility::vector1< Size > const & cutpoint_open = const_full_model_info( pose ).cutpoint_open_in_full_model();
		utility::vector1< Size > const & input_domain_map = const_full_model_info( pose ).input_domain_map();
		utility::vector1< Size > const & working_res = const_full_model_info( pose ).working_res();
		std::string const & full_sequence = const_full_model_info( pose ).full_sequence();
		missing_residues.clear();

		Size nmissing( 0 );
		Size const nres = pose_domain_map.size();
		for ( Size n = 1; n <= nres; n++ ){
			if ( input_domain_map[ n ] == 0 ){
				// missing residues.
				if ( pose_domain_map[ n ] == 0 && working_res.has_value( n ) ) {
					nmissing++;
					missing_residues.push_back( full_sequence[ n - 1 ] );
					}
			} else if ( n < nres &&
									!cutpoint_open.has_value( n ) &&
									input_domain_map[ n   ] > 0 && /* put in later, during loop_graph checks --rd2014 */
									input_domain_map[ n+1 ] > 0 &&
									input_domain_map[ n+1 ] != input_domain_map[ n ] ) {
				// missing suites between fixed domains (happens in, e.g., four-way junctions with no internal residues)
				if ( pose_domain_map[ n ] == 0 || /* don't need this? check with assert below */
						 pose_domain_map[ n+1 ] == 0 || /* don't need this? check with assert below */
						 ( pose_domain_map[ n ] != pose_domain_map[ n+1 ] ) ) {
					runtime_assert( pose_domain_map[ n ] != 0 );
					runtime_assert( pose_domain_map[ n+1 ] != 0 );
					nmissing++;
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////////
		// Testing use of LoopGraph. -- rhiju, dec. 2014
		//
		//  if runtime_assert()'s are not triggered, then deprecate above... in fact
		//   just carve out separate functions for get_number_missing_residues_and_connections,
		//   which is kind of a silly mixed up function.
		//
		//////////////////////////////////////////////////////////////////////////////
		core::scoring::loop_graph::LoopGraph loop_graph;
		loop_graph.update_loops( pose );
		runtime_assert( nmissing == loop_graph.nmissing( pose ) );
		runtime_assert( missing_residues == loop_graph.missing_residues( pose ) );  // if OK, switch to this.
		loop_suites = loop_graph.loop_suites( false /*include_free_loops*/ );

		return nmissing;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	bool
	check_all_residues_sampled( pose::Pose const & pose ){
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & input_domain_map =  full_model_info.input_domain_map();
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		for ( Size i = 1; i <= res_list.size(); i++ ) {
			if ( input_domain_map[ res_list[i] ] != 0 ) return false;
		}
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	std::map< std::pair< Size, Size >, std::pair< Size, Size > >
	get_preferred_jump_pair_for_docking_domains( FullModelInfo const & full_model_info ) {

		utility::vector1< Size > const dock_domain_map = full_model_info.dock_domain_map();

		// look for user-defined jumps between dock domains.
		// this does not actually have to be recalculated as its const for a full_model_parameters.
		std::map< std::pair< Size, Size >, std::pair< Size, Size > > preferred_jump_pair;
		utility::vector1< std::pair< Size, Size > > jump_pairs = full_model_info.full_model_parameters()->get_res_list_as_pairs( JUMP );
		for ( Size n = 1; n <= jump_pairs.size(); n++ ){
			Size const i_full = jump_pairs[ n ].first;
			Size const j_full = jump_pairs[ n ].second;
			runtime_assert( i_full < j_full );
			Size const dock_domain_i = dock_domain_map[ i_full ];
			Size const dock_domain_j = dock_domain_map[ j_full ];
			if ( dock_domain_i == dock_domain_j ) continue;

			// only one jump_res allowed for each dock_domain.
			runtime_assert( preferred_jump_pair.find( std::make_pair( dock_domain_i, dock_domain_j ) ) ==
											preferred_jump_pair.end() );
			preferred_jump_pair[ std::make_pair( dock_domain_i, dock_domain_j ) ] = std::make_pair( i_full, j_full );

			runtime_assert( preferred_jump_pair.find( std::make_pair( dock_domain_j, dock_domain_i ) ) ==
											preferred_jump_pair.end() );
			preferred_jump_pair[ std::make_pair( dock_domain_j, dock_domain_i ) ] = std::make_pair( j_full, i_full );
		}
		return preferred_jump_pair;
	}


	////////////////////////////////////////////////////////////////////////
	// go through pose and other_poses and figure out which chains
	// are connected in those poses.
	utility::vector1< std::pair< Size, Size > >
	get_chain_connections( pose::Pose const & pose ) {
		utility::vector1< Size > const chains_full = get_chains_full( pose );
		utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map_const( pose );
		utility::vector1< std::pair< Size, Size > > chain_connections;
		for ( Size n = 1; n <= max( pose_domain_map ); n++ ) {
			std::set< Size > chains_in_pose;
			for ( Size k = 1; k <= pose_domain_map.size(); k++ ) {
				if ( pose_domain_map[ k ] == n ) chains_in_pose.insert( chains_full[ k ] );
			}
			for ( std::set< Size >::const_iterator it1 = chains_in_pose.begin(), end = chains_in_pose.end(); it1 != end; ++it1 ) {
				for ( std::set< Size >::const_iterator it2 = it1; it2 != end; ++it2 ) {
					if ( it1 != it2 ) chain_connections.push_back( std::make_pair( *it1, *it2 ) );
				}
			}
		}
		return chain_connections;
	}

	/////////////////////////////////////////////////////////////////////////////////
	// take edges on graph and find connected clusters. I.e., the closure problem.
	// I think there's a faster way to do this, of course, and it probably exists elsewhere in the code,
	// but this doesn't need to be too fast.
	// yea, its in boost. if you're reading this, feel free to update.
	utility::vector1< Size >
	get_connection_domains( utility::vector1< std::pair< Size, Size > > const & chain_connections, Size const nchains ) {
		utility::vector1< Size > connection_domains;
		for ( Size k = 1; k <= nchains; k++ ) connection_domains.push_back( k );
		for ( Size n = 1; n <= chain_connections.size(); n++ ) {
			// coalesce domains that are connected by link.
			Size domain_old = std::max( connection_domains[ chain_connections[n].first ],
																	connection_domains[ chain_connections[n].second ] );
			Size domain_new = std::min( connection_domains[ chain_connections[n].first ],
																	connection_domains[ chain_connections[n].second ] );
			for ( Size k = 1; k <= nchains; k++ ) {
				if ( connection_domains[ k ]  == domain_old ) connection_domains[ k ] = domain_new;
			}
		}

		// need to relabel so that domains have numbers 1, 2, 3...
		Size count( 0 );
		std::set< Size > unique_domains( connection_domains.begin(), connection_domains.end() );
		utility::vector1< Size > connection_domains_relabel = connection_domains;
		for ( std::set< Size >::const_iterator it = unique_domains.begin(), end = unique_domains.end(); it != end; ++it ) {
			count++;
			Size const & domain_number( *it );
			for ( Size k = 1; k <= nchains; k++ ) {
				if ( connection_domains[ k ] == domain_number ) connection_domains_relabel[ k ] = count;
			}
		}

		return connection_domains_relabel;
	}

//////////////////////////////////////////////////////////////////////////////
bool
check_sample_sugar_in_full_model_info( pose::Pose const & pose,
																			 Size const i ) {
	if ( !full_model_info_defined( pose ) ) return false;
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & sample_sugar_res = full_model_info.rna_sample_sugar_res();
	return sample_sugar_res.has_value( full_model_info.res_list()[ i ] );
}


} //full_model_info
} //pose
} //core
