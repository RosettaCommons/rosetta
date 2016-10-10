// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/methods/carbohydrates/SugarBackboneEnergyCreator.hh
/// @brief   Method declarations for SugarBackboneEnergyCreator.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergyCreator_HH
#define INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergyCreator_HH

// Package headers
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>


namespace core {
namespace scoring {
namespace methods {
namespace carbohydrates {

/// @brief  EnergyMethodCreator allowing the ScoringManager to create a SugarBackboneEnergy method class
class SugarBackboneEnergyCreator : public EnergyMethodCreator {
public:
	/// @brief  Return an up-casted owning pointer (EnergyMethodOP) to the energy method.
	virtual EnergyMethodOP create_energy_method( EnergyMethodOptions const & ) const;

	/// @brief  Return the set of ScoreTypes for which this EnergyMethod is responsible.
	virtual ScoreTypes score_types_for_method() const;
};

}  // namespace carbohydrates
}  // namespace methods
}  // namespace scoring
}  // namespace core

#endif  // INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergyCreator_HH
