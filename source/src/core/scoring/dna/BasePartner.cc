// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <core/scoring/dna/BasePartner.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace dna {

/// helper routine

BasePartner const &
retrieve_base_partner_from_pose( pose::Pose const & pose )
{
	//using core::pose::datacache::CacheableDataType::BASE_PARTNER;

	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::BASE_PARTNER ) );
	debug_assert( dynamic_cast< BasePartner const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER ))));
	return ( static_cast< BasePartner const &>(    pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER )));
}


} // namespace dna
} // scoring
} // core

