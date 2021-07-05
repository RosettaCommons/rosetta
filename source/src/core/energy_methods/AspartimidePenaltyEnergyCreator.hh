// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/AspartimidePenaltyEnergyCreator.hh
/// @brief  Declaration for the class that connects AspartimidePenaltyEnergy with the ScoringManager.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_energy_methods_AspartimidePenaltyEnergyCreator_hh
#define INCLUDED_core_energy_methods_AspartimidePenaltyEnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>



namespace core {
namespace energy_methods {


class AspartimidePenaltyEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:
	/// @brief Instantiate a new AspartimidePenaltyEnergy
	core::scoring::methods::EnergyMethodOP
	create_energy_method(
		core::scoring::methods::EnergyMethodOptions const &options
	) const override;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	core::scoring::ScoreTypes
	score_types_for_method() const override;

};

}
}

#endif
