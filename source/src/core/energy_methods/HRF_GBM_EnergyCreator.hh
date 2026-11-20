// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/energy_methods/HRF_GBM_EnergyCreator.hh
/// @brief  Declaration of class that connects HRF_GBM_Energy with the ScoringManager
/// @author Elijah Day (ehday@ucla.edu)

#ifndef INCLUDED_core_energy_methods_HRF_GBM_EnergyCreator_hh
#define INCLUDED_core_energy_methods_HRF_GBM_EnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

namespace core {
namespace energy_methods {


class HRF_GBM_EnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:
	/// @brief Instantiate a new HRF_GBM_Energy
	core::scoring::methods::EnergyMethodOP
	create_energy_method(
		core::scoring::methods::EnergyMethodOptions const &
	) const override;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	core::scoring::ScoreTypes
	score_types_for_method() const override;

};

}
}

#endif // INCLUDED_core_energy_methods_HRF_GBM_EnergyCreator_hh
