// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/FaMPSolvEnergyCreator.hh
///
/// @brief  LK-Type Membrane Solvation Energy
/// @details Last Modified: 5/14/14
///
/// @author  Patrick Barth (Original)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_FaMPSolvEnergyCreator_hh
#define INCLUDED_core_scoring_membrane_FaMPSolvEnergyCreator_hh

// Unit Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Creator Class for Membrane CBeta Energy Method
class FaMPSolvEnergyCreator : public core::scoring::methods::EnergyMethodCreator {

public:

	/// @brief Instantiate a new MPEnvEnergy
	virtual
	methods::EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Return MPEnv Score Type Claimed by this energy method
	virtual
	ScoreTypes
	score_types_for_method() const;
};


} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_FaMPSolvEnergyCreator_hh
