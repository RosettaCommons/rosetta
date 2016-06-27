// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPNonHelixPenaltyCreator.hh
///
/// @brief  Membrane Protein Non helix in Mmebrane Penalty
/// @details Whole structure energy - penalty for helices not in the membrane?
///    and uses mpframework data
///    Last Modified: 3/31/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPNonHelixPenaltyCreator_hh
#define INCLUDED_core_scoring_membrane_MPNonHelixPenaltyCreator_hh

// Unit Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Membrane Non Hleix in mmebrane Penalty Creator Class
class MPNonHelixPenaltyCreator : public methods::EnergyMethodCreator
{
public:

	/// @brief Instantiate a new MPNonHelixPenalty
	virtual
	methods::EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	ScoreTypes
	score_types_for_method() const;

};

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MPNonHelixPenaltyCreator_hh
