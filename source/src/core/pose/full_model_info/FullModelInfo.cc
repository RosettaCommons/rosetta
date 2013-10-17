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
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>

// Package headers
#include <core/pose/full_model_info/FullModelInfoUtil.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/vector1.hh>

// C++
#include <string>
#include <map>

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
FullModelInfo::FullModelInfo() {} // blank. Should not be used?

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FullModelInfo::FullModelInfo( std::string const full_sequence ):
	CacheableData(),
	full_sequence_( full_sequence )
{
	for ( Size n = 1; n <= full_sequence.size(); n++ ) 		conventional_numbering_.push_back( n );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FullModelInfo::FullModelInfo( pose::Pose & pose,
															std::string const & full_sequence,
															utility::vector1< Size > const & cutpoint_open_in_full_model,
															utility::vector1< Size > const & res_numbers_in_pose ):
	CacheableData(),
	full_sequence_( full_sequence ),
	cutpoint_open_in_full_model_( cutpoint_open_in_full_model ),
	res_list_( res_numbers_in_pose )
{
	for ( Size n = 1; n <= full_sequence.size(); n++ ){
		if ( res_numbers_in_pose.has( n ) ){
			fixed_domain_map_.push_back( 1 );
		} else {
			fixed_domain_map_.push_back( 0 );
		}
		conventional_numbering_.push_back( n );
	}
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FullModelInfo::FullModelInfo( pose::Pose & pose ) :
	CacheableData()
{

	res_list_ = get_res_num_from_pdb_info( pose );

	get_sequence_with_gaps_filled_with_n( pose, full_sequence_, conventional_numbering_ );
	cutpoint_open_in_full_model_ = get_cutpoint_open_from_pdb_info( pose );

	// not sure what's best here -- for now setting that the pose's residues are 'fixed' within the domain map.
	fixed_domain_map_ = utility::vector1<Size>( full_sequence_.size(), 0 );

	for ( Size n = 1; n <= res_list_.size(); n++ ) {
		Size const & res_num = res_list_[ n ];
		runtime_assert( conventional_numbering_.has_value( res_num ) );
		fixed_domain_map_[ conventional_numbering_.index( res_num ) ] = 1;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Copy constructors must copy all data, not just some...
FullModelInfo::FullModelInfo( FullModelInfo const & src ) :
	CacheableData(),
	full_sequence_( src.full_sequence_ ),
	cutpoint_open_in_full_model_( src.cutpoint_open_in_full_model_ ),
	fixed_domain_map_( src.fixed_domain_map_ ),
	conventional_numbering_( src.conventional_numbering_ ),
	res_list_( src.res_list_ )
{
	// we have to tell our daughters in the pose tree that its time to get cloned.
	for ( Size n = 1; n <= src.other_pose_list_.size(); n++ ){
		other_pose_list_.push_back( src.other_pose_list_[ n ]->clone() );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FullModelInfo::~FullModelInfo(){}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::get_res_num_from_pdb_info( pose::Pose const & pose ) const {

	utility::vector1< Size > resnum;

	PDBInfoCOP pdb_info = pose.pdb_info();

	if ( pdb_info )	{
		for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( pdb_info->number( n ) );
	} else {
		for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( n );
	}

	return resnum;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::get_sequence_with_gaps_filled_with_n( pose::Pose const & pose,
																										 std::string & sequence,
																										 utility::vector1< Size > & full_numbering ) const {

	// should also be smart about not filling in n's between chains.
	// anyway. this is a quick hack for now.
	utility::vector1< Size > const & res_list = res_list_;

	sequence = "";
	full_numbering.clear();

	sequence.push_back( pose.sequence()[ 0 ]  );
	full_numbering.push_back( res_list[ 1 ] );

	for ( Size n = 2; n <= pose.total_residue(); n++ ){

		Size const prev_res_num    = res_list[ n-1 ];
		Size const current_res_num = res_list[ n ];
		for ( Size i = prev_res_num+1; i < current_res_num; i++ ) {
			sequence.push_back( 'n' );
			full_numbering.push_back( i );
		}
		sequence.push_back( pose.sequence()[ n-1 ] );
		full_numbering.push_back( current_res_num );
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::chains_in_full_model() const {
	utility::vector1< Size > chains;
	Size chain_number( 1 );
	for ( Size n = 1; n <= full_sequence_.size(); n++ ){
		chains.push_back( chain_number );
		if ( cutpoint_open_in_full_model_.has_value( n ) ) chain_number++;
	}
	return chains;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::moving_res_in_full_model() const {
	utility::vector1< Size > moving_res;
	runtime_assert( full_sequence_.size() == fixed_domain_map_.size() );
	for ( Size n = 1; n <= full_sequence_.size(); n++ ){
		if ( fixed_domain_map_[n] == 0 ) moving_res.push_back( n );
	}
	return moving_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const {

	PDBInfoCOP pdb_info = pose.pdb_info();
	utility::vector1< Size > cutpoint_open;

	for ( Size n = 1; n < pose.total_residue(); n++ ){

		if ( pdb_info &&  (pdb_info->chain( n ) != pdb_info->chain( n+1 )) )	{
			cutpoint_open.push_back( n );
			continue;
		}

		//		if ( pose.fold_tree().is_cutpoint( n ) &&
		//				 !pose.residue( n ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
		//				 !pose.residue( n+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
		//			cutpoint_open.push_back( n );
		//		}
	}

	return cutpoint_open;
}

///////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::full_to_sub( utility::vector1< Size > const & res_in_full_model_numbering ) const{
	utility::vector1< Size > res;
	for ( Size n = 1; n <= res_in_full_model_numbering.size(); n++ ){
		res.push_back( res_list_.index( res_in_full_model_numbering[ n ] ) );
	}
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::find_index_in_other_pose_list( pose::Pose const & pose) const {

	for ( Size n = 1; n <= other_pose_list_.size(); n++ ){
		if ( other_pose_list_[ n ] == & pose ) return n;
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::clear_other_pose_list() {
	other_pose_list_.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::get_idx_for_other_pose_with_residue( Size const input_res ) const {
	for ( Size i = 1; i <= other_pose_list_.size(); i++ ){
		utility::vector1< Size > const & daughter_res_list = const_full_model_info( *other_pose_list_[i] ).res_list();
		if ( daughter_res_list.has_value( input_res ) ) return i;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::add_other_pose( core::pose::PoseOP & pose )
{
	other_pose_list_.push_back( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::set_other_pose_list( utility::vector1< pose::PoseOP > const & setting ){
	other_pose_list_ = setting;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::remove_other_pose_at_idx( Size const idx ){

	runtime_assert( idx <= other_pose_list_.size() );
	utility::vector1< core::pose::PoseOP > other_pose_list_new;

	for ( Size i = 1; i <= other_pose_list_.size(); i++ ) {
		if ( i == idx ) continue;
		other_pose_list_new.push_back( other_pose_list_[ i ] );
	}

	other_pose_list_ = other_pose_list_new;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// @details Pose must already contain a full_model_info object or this method will fail.
FullModelInfo const &
const_full_model_info( pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) );
	return *( static_cast< FullModelInfo const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO)() ) );
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Either returns a non-const reference to the FullModelInfo object already stored
/// in the pose, or creates a new FullModelInfo object, places it in the pose, and returns
/// a non-const reference to it.
FullModelInfo &
nonconst_full_model_info( pose::Pose & pose )
{

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) ) {
		// following takes time -- so we shouldn't call this nonconst function a lot.
		//		if ( check_full_model_info_OK( pose ) ){
		runtime_assert ( check_full_model_info_OK( pose ) ); // later remove this for speed.
		return *( static_cast< FullModelInfo * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO )() ));
			//		}
	}

	FullModelInfoOP full_model_info = new FullModelInfo( pose );

	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
	return *full_model_info;
}

}
}
}
