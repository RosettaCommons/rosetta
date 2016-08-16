// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/methods/carbohydrates/SugarBackboneEnergy.fwd.hh
/// @brief   Forward declarations for SugarBackboneEnergy.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergy_FWD_HH
#define INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergy_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace methods {
namespace carbohydrates {

/// @brief  An energy method class for scoring a carbohydrate backbone based on its phi and psi angles.
class SugarBackboneEnergy;

typedef utility::pointer::shared_ptr< SugarBackboneEnergy > SugarBackboneEnergyOP;
typedef utility::pointer::shared_ptr< SugarBackboneEnergy const> SugarBackboneEnergyCOP;

}  // namespace carbohydrates
}  // namespace methods
}  // namespace scoring
}  // namespace core

#endif // INCLUDED_core_scoring_methods_carbohydrates_SugarBackboneEnergy_FWD_HH
