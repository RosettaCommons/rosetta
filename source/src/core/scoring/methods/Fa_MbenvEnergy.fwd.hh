// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/methods/Fa_MbenvEnergy.fwd.hh
/// @brief  Probability of observing an amino acid (NOT conditional on phi/psi), energy method forward declaration
/// @author Ron Jacak


#ifndef INCLUDED_core_scoring_methods_Fa_MbenvEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_Fa_MbenvEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class Fa_MbenvEnergy;

typedef utility::pointer::shared_ptr< Fa_MbenvEnergy > Fa_MbenvEnergyOP;

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_methods_Fa_MbenvEnergy_FWD_HH
