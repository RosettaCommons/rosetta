// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/SplitUnfoldedTwoBodyEnergyCreator.hh
/// @brief  Energy creator for the split unfolded two body energy method
/// @author Riley Simmons-Edler (rse231@nyu.edu)

#ifndef INCLUDED_core_energy_methods_SplitUnfoldedTwoBodyEnergyCreator_hh
#define INCLUDED_core_energy_methods_SplitUnfoldedTwoBodyEnergyCreator_hh

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>




namespace core {
namespace energy_methods {


class SplitUnfoldedTwoBodyEnergyCreator : public core::scoring::methods::EnergyMethodCreator
{
public:
	core::scoring::methods::EnergyMethodOP create_energy_method(const core::scoring::methods::EnergyMethodOptions &) const override;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	core::scoring::ScoreTypes score_types_for_method() const override;
};

}
}

#endif
