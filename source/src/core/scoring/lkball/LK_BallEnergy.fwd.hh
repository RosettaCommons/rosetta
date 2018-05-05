// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LK_BallEnergy.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author David Baker
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_LK_BALLENERGY_FWD_HH
#define INCLUDED_core_scoring_methods_LK_BALLENERGY_FWD_HH

// Unit Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace lkball {

class LK_BallEnergy;

typedef utility::pointer::shared_ptr< LK_BallEnergy > LK_BallEnergyOP;
typedef utility::pointer::shared_ptr< LK_BallEnergy const > LK_BallEnergyCOP;

class LK_Ball_RPE_Invoker;

}
}
}

#endif // INCLUDED_core_scoring_methods_LK_BallEnergy_HH
