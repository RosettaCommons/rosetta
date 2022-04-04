// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/lkball/LK_DomeEnergy.fwd.hh
/// @brief  lk_dome second shell water interactions
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_scoring_lkball_LK_DomeEnergy_fwd_hh
#define INCLUDED_core_scoring_lkball_LK_DomeEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace lkball {


class LK_DomeEnergy;

typedef utility::pointer::shared_ptr< LK_DomeEnergy > LK_DomeEnergyOP;
typedef utility::pointer::shared_ptr< LK_DomeEnergy const > LK_DomeEnergyCOP;


} // methods
} // scoring
} // core

#endif
