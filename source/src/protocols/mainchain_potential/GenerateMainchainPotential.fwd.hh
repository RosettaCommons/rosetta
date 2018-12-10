// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mainchain_potential/GenerateMainchainPotential.fwd.hh
/// @brief Forward declarations for a generator for mainchain potentials.  Inputs are a noncanonical residue type with an already-generated
/// sidechain potential; outputs are a potential file suitable for use by the RamaPrePro scoreterm.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_mainchain_potential_GenerateMainchainPotential_fwd_hh
#define INCLUDED_protocols_mainchain_potential_GenerateMainchainPotential_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace mainchain_potential {

class GenerateMainchainPotential;

typedef utility::pointer::shared_ptr< GenerateMainchainPotential > GenerateMainchainPotentialOP;
typedef utility::pointer::shared_ptr< GenerateMainchainPotential const > GenerateMainchainPotentialCOP;

} //protocols
} //mainchain_potential

#endif //INCLUDED_protocols_mainchain_potential_GenerateMainchainPotential_fwd_hh
