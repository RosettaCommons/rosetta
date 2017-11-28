// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/netcharge_energy/NetChargeEnergySetup.fwd.hh
/// @brief Forward declarations for a helper class that stores the setup information for the NetChargeEnergy score term.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_core_scoring_netcharge_energy_NetChargeEnergySetup_fwd_hh
#define INCLUDED_core_scoring_netcharge_energy_NetChargeEnergySetup_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace netcharge_energy {

class NetChargeEnergySetup;

typedef utility::pointer::shared_ptr< NetChargeEnergySetup > NetChargeEnergySetupOP;
typedef utility::pointer::shared_ptr< NetChargeEnergySetup const > NetChargeEnergySetupCOP;

} // netcharge_energy
} // scoring
} // core


#endif
