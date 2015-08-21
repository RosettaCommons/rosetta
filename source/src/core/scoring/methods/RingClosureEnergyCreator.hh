// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RingClosureEnergyCreator.hh
/// @brief  Declaration for the class that connects RingClosureEnergy with the ScoringManager
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory

#ifndef INCLUDED_core_scoring_methods_RingClosureEnergyCreator_hh
#define INCLUDED_core_scoring_methods_RingClosureEnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class RingClosureEnergyCreator : public EnergyMethodCreator
{
public:
	/// @brief Instantiate a new RingClosureEnergy
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

} //methods
} //scoring
} //core

#endif //INCLUDED_core_scoring_methods_RingClosureEnergyCreator_hh
