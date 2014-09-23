// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SA_Energy.cc
/// @brief  sa energy function definition.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/SA_Energy.hh>
#include <core/scoring/methods/SA_EnergyCreator.hh>

// Package headers
#include <core/scoring/sasa.hh>
// AUTO-REMOVED #include <core/scoring/packing/surf_vol.hh>
// Project headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the SA_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
SA_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new SA_Energy );
}

ScoreTypes
SA_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( sa );
	return sts;
}

/// c-tor
SA_Energy::SA_Energy() :
	parent( methods::EnergyMethodCreatorOP( new SA_EnergyCreator ) )
{}


/// clone
EnergyMethodOP
SA_Energy::clone() const
{
	return EnergyMethodOP( new SA_Energy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
SA_Energy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace conformation;

	totals[ sa ] = calc_total_sasa( pose, 1.4 ); //default water probe
//	totals[ sa ] = core::scoring::packing::get_surf_tot(pose, 1.4); //default water probe

} // finalize_total_energy
core::Size
SA_Energy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core
