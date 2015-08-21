// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/AARepeatEnergyCreator.hh
/// @brief Creator for an EnergyMethod that penalizes stretches of a repeating amino acid (e.g. poly-Q sequences).
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_scoring_methods_AARepeatEnergyCreator_hh
#define INCLUDED_core_scoring_methods_AARepeatEnergyCreator_hh

// Unit header
#include <core/scoring/methods/EnergyMethodCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class AARepeatEnergyCreator : public EnergyMethodCreator
{
public:

	/// @brief Instantiate a new AARepeatEnergy.
	///
	virtual methods::EnergyMethodOP create_energy_method( methods::EnergyMethodOptions const & ) const;

	/// @brief Return the set of score types claimed by the EnergyMethod that
	/// this EnergyMethodCreator creates in its create_energy_method() function.
	virtual ScoreTypes score_types_for_method() const;
};

} // methods
} // scoring
} // core

#endif
