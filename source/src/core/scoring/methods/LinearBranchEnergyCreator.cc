// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LinearBranchEnergyCreator.cc
/// @brief  LinearBranchEnergyCreator implementation
/// @author Christopher Miles (cmiles@uw.edu)

#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/LinearBranchEnergy.hh>
#include <core/scoring/methods/LinearBranchEnergyCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @details Creates and initializes a new LinearBranchEnergy object with the
/// specified options
EnergyMethodOP
LinearBranchEnergyCreator::create_energy_method(const EnergyMethodOptions& /*opt*/) const {
	return EnergyMethodOP( new LinearBranchEnergy );
}

ScoreTypes
LinearBranchEnergyCreator::score_types_for_method() const {
	ScoreTypes types;
	types.push_back(linear_branch_conn);
	return types;
}

}  // namespace methods
}  // namespace scoring
}  // namespace core
