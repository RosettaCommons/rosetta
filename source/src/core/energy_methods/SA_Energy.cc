// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/SA_Energy.cc
/// @brief  sa energy function definition.
/// @author James Thompson


// Unit headers
#include <core/energy_methods/SA_Energy.hh>
#include <core/energy_methods/SA_EnergyCreator.hh>

// Package headers
#include <core/scoring/sasa.hh>
// Project headers

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the SA_Energy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
SA_EnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< SA_Energy >();
}

core::scoring::ScoreTypes
SA_EnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( sa );
	return sts;
}

/// c-tor
SA_Energy::SA_Energy() :
	parent( utility::pointer::make_shared< SA_EnergyCreator >() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
SA_Energy::clone() const
{
	return utility::pointer::make_shared< SA_Energy >();
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
SA_Energy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const {
	using namespace conformation;

	totals[ core::scoring::sa ] = core::scoring::calc_total_sasa( pose, 1.4 ); //default water probe

} // finalize_total_energy
core::Size
SA_Energy::version() const
{
	return 1; // Initial versioning
}

} // scoring
} // core
