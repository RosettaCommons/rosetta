// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/fiber_diffraction/FAScatter.cc
/// @brief Cache for scattering factors in all-atom mode
/// @author Wojciech Potrzebowski and Ingemar Andre

#include <core/scoring/fiber_diffraction/FAScatter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


namespace core {
namespace scoring {
namespace fiber_diffraction {

/// helper routine

FAScatter &
retrieve_fa_scatter_from_pose( pose::Pose & pose )
{
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_FA_SCATTERING ) );
	debug_assert( dynamic_cast< FAScatter *>( &( pose.data().get( core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_FA_SCATTERING ))));
	return ( static_cast< FAScatter &>(    pose.data().get( core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_FA_SCATTERING )));
	//return utility::pointer::static_pointer_cast< CentroidScatter > ( pose.data().get_ptr(
	//      core::pose::datacache::CacheableDataType::FIBER_DIFFRACTION_CEN_SCATTERING) );

}


} // namespace fiber_diffraction
} // scoring
} // core

