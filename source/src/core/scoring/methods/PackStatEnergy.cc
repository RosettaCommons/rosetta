// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PackStatEnergy.cc
/// @brief  Radius of gyration energy function definition. Returns -1 * RG for a given Pose.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/PackStatEnergy.hh>
#include <core/scoring/methods/PackStatEnergyCreator.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/packstat/compute_sasa.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the PackStatEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
PackStatEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new PackStatEnergy;
}

ScoreTypes
PackStatEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( pack_stat );
	return sts;
}


/// c-tor
PackStatEnergy::PackStatEnergy() :
	parent( methods::EnergyMethodCreatorOP( new PackStatEnergyCreator ) )
{}


/// clone
EnergyMethodOP
PackStatEnergy::clone() const
{
	return new PackStatEnergy();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the radius of gyration and place the answer into
/// totals[ rg ].
void
PackStatEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace core::scoring::packstat;

	core::Real packing_score = compute_packing_score(pose);
	// std::cerr << "ps " << packing_score << " " << std::endl;
	totals[ pack_stat ] = - packing_score * pose.total_residue();

}
core::Size
PackStatEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
