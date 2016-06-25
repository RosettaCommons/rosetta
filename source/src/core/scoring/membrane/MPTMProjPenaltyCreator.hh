// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPTMProjPenaltyCreator.hh
///
/// @brief  Membrane Protein TM Proj Penalty
/// @details Whole structure energy - Penalty for unreasonable tm-helix length compared to predicted
///    helix length (from topology) and uses mpframework data
///    Last Modified: 4/3/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPTMProjPenaltyCreator_hh
#define INCLUDED_core_scoring_membrane_MPTMProjPenaltyCreator_hh

// Unit Headers
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>

// Rebecca, our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::scoring;

namespace core {
namespace scoring {
namespace membrane {

/// @brief Membrane TM Proj Penalty Creator Class
class MPTMProjPenaltyCreator : public methods::EnergyMethodCreator
{
public:

	/// @brief Instantiate a new MPTMProjPenalty
	virtual
	methods::EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Returns MPTMProj score type
	virtual
	ScoreTypes
	score_types_for_method() const;

};

} // membrane
} // scoring
} // core


#endif // INCLUDED_core_scoring_membrane_MPTMProjPenaltyCreator_hh
