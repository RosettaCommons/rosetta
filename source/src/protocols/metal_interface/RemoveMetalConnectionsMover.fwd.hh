// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/metal_interface/RemoveMetalConnectionsMover.fwd.hh
/// @brief A mover that removes the connections to metals that were added by the SetupMetalsMover
/// or by the -auto_setup_metals flag.
/// @details This mover:
///     - Removes the bonds between metals and metal-binding residues.
///     - Reverts metal-liganding residues back to their pre-bonded types.
///     - Reverts metals back to their pre-bonded types.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMover_fwd_hh
#define INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace metal_interface {

class RemoveMetalConnectionsMover;

using RemoveMetalConnectionsMoverOP = utility::pointer::shared_ptr< RemoveMetalConnectionsMover >;
using RemoveMetalConnectionsMoverCOP = utility::pointer::shared_ptr< RemoveMetalConnectionsMover const >;

} //metal_interface
} //protocols

#endif //INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMover_fwd_hh
