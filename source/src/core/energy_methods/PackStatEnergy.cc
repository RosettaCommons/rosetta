// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/PackStatEnergy.cc
/// @brief  Radius of gyration energy function definition. Returns -1 * RG for a given Pose.
/// @author James Thompson


// Unit headers
#include <core/energy_methods/PackStatEnergy.hh>
#include <core/energy_methods/PackStatEnergyCreator.hh>

// Package headers
#include <core/scoring/packstat/compute_sasa.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the PackStatEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
PackStatEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< PackStatEnergy >();
}

core::scoring::ScoreTypes
PackStatEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( pack_stat );
	return sts;
}

/// c-tor
PackStatEnergy::PackStatEnergy() :
	parent( utility::pointer::make_shared< PackStatEnergyCreator >() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
PackStatEnergy::clone() const
{
	return utility::pointer::make_shared< PackStatEnergy >();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the radius of gyration and place the answer into
/// totals[ core::scoring::rg ].
void
PackStatEnergy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const {
	using namespace core::scoring::packstat;

	core::Real packing_score = compute_packing_score(pose);
	// std::cerr << "ps " << packing_score << " " << std::endl;
	totals[ core::scoring::pack_stat ] = - packing_score * pose.size();
}

core::Size
PackStatEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
