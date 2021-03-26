// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaProtocol_v1.fwd.hh
/// @brief Version 1 model for trRosetta.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_trRosettaProtocol_v1_fwd_hh
#define INCLUDED_protocols_trRosetta_trRosettaProtocol_v1_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace trRosetta {

class trRosettaProtocol_v1;

using trRosettaProtocol_v1OP = utility::pointer::shared_ptr< trRosettaProtocol_v1 >;
using trRosettaProtocol_v1COP = utility::pointer::shared_ptr< trRosettaProtocol_v1 const >;

} //trRosetta
} //protocols

#endif //INCLUDED_protocols_trRosetta_trRosettaProtocol_v1_fwd_hh
