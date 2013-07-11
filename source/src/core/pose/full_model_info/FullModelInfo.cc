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
FullModelInfo::FullModelInfo( 	utility::vector1< Size > const & sub_to_full,
																utility::vector1< Size > const & moving_res_list,
																std::string const full_sequence,
																utility::vector1< Size > const & cutpoint_open_in_full_model) : //proper constructor
	CacheableData(),
	sub_to_full_( sub_to_full ),
	moving_res_list_( moving_res_list ),
	full_sequence_( full_sequence ),
	cutpoint_open_in_full_model_( cutpoint_open_in_full_model )
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FullModelInfo::FullModelInfo( pose::Pose const & pose ) :
	CacheableData(),
	sub_to_full_( get_res_num_from_pdb_info( pose ) ),
	full_sequence_( pose.sequence() ),
	cutpoint_open_in_full_model_( get_cutpoint_open_from_pdb_info( pose ) )
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Copy constructors must copy all data, not just some...
FullModelInfo::FullModelInfo( FullModelInfo const & src ) :
	CacheableData(),
	sub_to_full_( src.sub_to_full_ ),
	moving_res_list_( src.moving_res_list_ ),
	full_sequence_( src.full_sequence_ ),
	cutpoint_open_in_full_model_( src.cutpoint_open_in_full_model_ )
{
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::get_res_num_from_pdb_info( pose::Pose const & pose ) const{

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
utility::vector1< Size >
FullModelInfo::get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const{

	PDBInfoCOP pdb_info = pose.pdb_info();
	utility::vector1< Size > cutpoint_open;

	if ( pdb_info )	 {

		for ( Size n = 1; n < pose.total_residue(); n++ ){
			if ( pdb_info->chain( n ) != pdb_info->chain( n+1 ) )	cutpoint_open.push_back( n ); // hopefully this will be smart about cutpoints? no it isn't.
		}

	} else {

		for ( Size n = 1; n < pose.total_residue(); n++ ){
			if ( pose.fold_tree().is_cutpoint( n ) &&
					 !pose.residue( n ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
					 !pose.residue( n+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
				cutpoint_open.push_back( n );
			}
		}

	}

	return cutpoint_open;
}

/// @details Pose must already contain a full_model_info object or this method will fail.
FullModelInfo const &
const_full_model_info_from_pose( pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) );
	return *( static_cast< FullModelInfo const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO)() ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Either returns a non-const reference to the rna_scoring object already stored
/// in the pose, or creates a new rna scoring info object, places it in the pose, and returns
/// a non-const reference to it.
FullModelInfo &
nonconst_full_model_info_from_pose( pose::Pose & pose )
{
	if ( pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) ) {
		return *( static_cast< FullModelInfo * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO )() ));
	}

	FullModelInfoOP full_model_info = new FullModelInfo( pose );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
	return *full_model_info;
}

}
}
}
