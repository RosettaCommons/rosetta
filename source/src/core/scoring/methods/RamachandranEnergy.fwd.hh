// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RamachandranEnergy.fwd.hh
/// @brief  Ramachandran energy method class forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_methods_RamachandranEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_RamachandranEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class RamachandranEnergy;

typedef utility::pointer::shared_ptr< RamachandranEnergy > RamachandranEnergyOP;
typedef utility::pointer::shared_ptr< RamachandranEnergy const > RamachandranEnergyCOP;

} // methods
} // scoring
} // core


#endif
