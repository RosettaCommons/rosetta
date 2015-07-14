// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RingClosureEnergy.fwd.hh
/// @brief  Noncanonical ring closure energy method class forward declaration
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory


#ifndef INCLUDED_core_scoring_methods_RingClosureEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_RingClosureEnergy_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class RingClosureEnergy;

typedef utility::pointer::shared_ptr< RingClosureEnergy > RingClosureEnergyOP;
typedef utility::pointer::shared_ptr< RingClosureEnergy const > RingClosureEnergyCOP;

} // methods
} // scoring
} // core

#endif
