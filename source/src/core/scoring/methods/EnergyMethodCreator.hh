// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EnergyMethodCreator.hh
/// @brief  Declaration of the base class for EnergyMethod factory registration and creation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_methods_EnergyMethodCreator_hh
#define INCLUDED_core_scoring_methods_EnergyMethodCreator_hh

// Unit headers
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief The EnergyMethodCreator class's responsibilities are to create
/// on demand a new EnergyMethod class, and to tell the ScoringManager
/// singleton which ScoreTypes the EnergyMethod it creates is responsible for.
/// The EnergyMethodCreator must register itself with the ScoringManager
/// at load time (before main() begins) so that the ScoringManager is ready
/// to start creating EnergyMethods by the time the first ScoreFunction
/// requests one.
class EnergyMethodCreator : public utility::pointer::ReferenceCount
{
public:
	/// @brief Instantiate a new EnergyMethod given a set of energy-method options
	virtual
	EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const & options
	) const = 0;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	ScoreTypes
	score_types_for_method() const = 0;

};

}
}
}

#endif
