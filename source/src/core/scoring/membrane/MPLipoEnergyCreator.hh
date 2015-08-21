// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/scoring/membrane/MPLipoEnergy.hh
///
/// @brief  Membrane Lipophibicity Term
/// @details Whole Structure Energy - Evaluate structure based on derived
///    lipophobicities from input in lips file.
///    Last Modified: 3/28/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MPLipoEnergyCreator_hh
#define INCLUDED_core_scoring_membrane_MPLipoEnergyCreator_hh

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

using namespace core::scoring;
using namespace core::scoring::methods;

/// @brief Creator Class for MPLipoEnergy
class MPLipoEnergyCreator : public EnergyMethodCreator
{
public:

	/// @brief Instantiate a new MPLipo
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

#endif // INCLUDED_core_scoring_membrane_MPLipoEnergyCreator_hh
