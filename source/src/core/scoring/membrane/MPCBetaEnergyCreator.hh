// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/membrane/MPCbetaEnergyCreator.hh
///
/// @brief  Membrane Environemnt CBeta Energy
/// @details One Body Term - Score packing density in the membrane. Scores centroids for within
///    6A and 12A radius. Derived from Membrane base potential and uses mpframework data
///    Last Modified: 4/2/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPCbetaEnergyCreator_hh
#define INCLUDED_core_scoring_membrane_MPCbetaEnergyCreator_hh

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
class MPCbetaEnergyCreator : public core::scoring::methods::EnergyMethodCreator {

public:

	/// @brief Instantiate a new MPCBetaEnergy
	virtual
	methods::EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Return the Set of Score Types claimed by EnergyMethod (uses MPCbeta)
	virtual
	ScoreTypes
	score_types_for_method() const;
};


} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPCbetaEnergyCreator_hh
