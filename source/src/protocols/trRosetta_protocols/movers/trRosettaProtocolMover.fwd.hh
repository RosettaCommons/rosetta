// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file trRosetta_protocols/movers/trRosettaProtocolMover.fwd.hh
/// @brief The full trRosetta structure prediction protocol from Yang et al, converted to
/// C++ and implemented as a mover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMover_fwd_hh
#define INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace trRosetta_protocols {
namespace movers {

class trRosettaProtocolMover;

using trRosettaProtocolMoverOP = utility::pointer::shared_ptr< trRosettaProtocolMover >;
using trRosettaProtocolMoverCOP = utility::pointer::shared_ptr< trRosettaProtocolMover const >;

} //movers
} //trRosetta_protocols
} //protocols

#endif //INCLUDED_protocols_trRosetta_protocols_movers_trRosettaProtocolMover_fwd_hh
