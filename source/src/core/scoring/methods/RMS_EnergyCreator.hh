// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/methods/RMS_EnergyCreator.hh
/// @brief  Declaration for the class that connects RMS_Energy with the ScoringManager
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_methods_RMS_EnergyCreator_hh
#define INCLUDED_core_scoring_methods_RMS_EnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class RMS_EnergyCreator : public EnergyMethodCreator
{
public:
	/// @brief Instantiate a new RMS_Energy
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

}
}
}

#endif
