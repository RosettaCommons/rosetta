// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mainchain_potential/GenerateMainchainPotentialOptions.fwd.hh
/// @brief Options container for the generator for mainchain potentials.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_mainchain_potential_GenerateMainchainPotentialOptions_fwd_hh
#define INCLUDED_protocols_mainchain_potential_GenerateMainchainPotentialOptions_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace mainchain_potential {

class GenerateMainchainPotentialOptions;

typedef utility::pointer::shared_ptr< GenerateMainchainPotentialOptions > GenerateMainchainPotentialOptionsOP;
typedef utility::pointer::shared_ptr< GenerateMainchainPotentialOptions const > GenerateMainchainPotentialOptionsCOP;

} //protocols
} //mainchain_potential

#endif //INCLUDED_protocols_mainchain_potential_GenerateMainchainPotentialOptions_fwd_hh
