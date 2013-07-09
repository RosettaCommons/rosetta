// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/swa/SubToFullInfo.cc
/// @brief  Mapping from a working pose into a bigger pose, for swa monte carlo stuff.
/// @author Rhiju Das

// Unit headers
#include <protocols/swa/monte_carlo/SubToFullInfo.hh>
#include <protocols/swa/monte_carlo/SubToFullInfo.fwd.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
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

namespace protocols {
namespace swa {
namespace monte_carlo {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SubToFullInfo::SubToFullInfo() {} // blank. Should not be used?

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SubToFullInfo::SubToFullInfo( 	std::map< Size, Size > sub_to_full,
																utility::vector1< Size > moving_res_list,
																std::string full_sequence,
																utility::vector1< Size > cutpoints_in_full_pose) : //proper constructor
	CacheableData(),
	sub_to_full_( sub_to_full ),
	moving_res_list_( moving_res_list ),
	full_sequence_( full_sequence ),
	cutpoints_in_full_pose_( cutpoints_in_full_pose )
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Copy constructors must copy all data, not just some...
	SubToFullInfo::SubToFullInfo( SubToFullInfo const & src ) :
	CacheableData(),
	sub_to_full_( src.sub_to_full_ ),
	moving_res_list_( src.moving_res_list_ ),
	full_sequence_( src.full_sequence_ ),
	cutpoints_in_full_pose_( src.cutpoints_in_full_pose_ )
{
}

/// @details Pose must already contain a sub_to_full_info object or this method will fail.
SubToFullInfo const &
const_sub_to_full_info_from_pose( pose::Pose & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::SUB_TO_FULL_INFO ) );
	return *( static_cast< SubToFullInfo const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::SUB_TO_FULL_INFO)() ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Either returns a non-const reference to the rna_scoring object already stored
/// in the pose, or creates a new rna scoring info object, places it in the pose, and returns
/// a non-const reference to it.
SubToFullInfo &
nonconst_sub_to_full_info_from_pose( pose::Pose & pose )
{
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SUB_TO_FULL_INFO ) ) {
		return *( static_cast< SubToFullInfo * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SUB_TO_FULL_INFO )() ));
	}

	SubToFullInfoOP sub_to_full_info = new SubToFullInfo();
	pose.data().set( core::pose::datacache::CacheableDataType::SUB_TO_FULL_INFO, sub_to_full_info );
	return *sub_to_full_info;
}

}
}
}
