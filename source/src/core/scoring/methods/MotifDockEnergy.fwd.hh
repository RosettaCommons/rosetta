// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MotifDockEnergy.hh
/// @brief  Adaptation of Motif score for Docking
/// @author Nick Marze (nickmarze@gmail.com)


#ifndef INCLUDED_core_scoring_methods_carbohydrates_MotifDockEnergy_FWD_HH
#define INCLUDED_core_scoring_methods_carbohydrates_MotifDockEnergy_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief  Adaptation of Motif score for Docking
class MotifDockEnergy;

typedef utility::pointer::shared_ptr< MotifDockEnergy > MotifDockEnergyOP;
typedef utility::pointer::shared_ptr< MotifDockEnergy const> MotifDockEnergyCOP;

}  // namespace methods
}  // namespace scoring
}  // namespace core

#endif // INCLUDED_core_scoring_methods_carbohydrates_MotifDockEnergy_FWD_HH
