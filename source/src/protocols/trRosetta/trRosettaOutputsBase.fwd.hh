// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaOutputsBase.fwd.hh
/// @brief A pure virtual base class for the outputs of trRosetta.  Derived classes are for particular trRosetta versions (to allow for future versions providing additional outputs).
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_trRosetta_trRosettaOutputsBase_fwd_hh
#define INCLUDED_protocols_trRosetta_trRosettaOutputsBase_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace trRosetta {

class trRosettaOutputsBase;

using trRosettaOutputsBaseOP = utility::pointer::shared_ptr< trRosettaOutputsBase >;
using trRosettaOutputsBaseCOP = utility::pointer::shared_ptr< trRosettaOutputsBase const >;

} //trRosetta
} //protocols

#endif //INCLUDED_protocols_trRosetta_trRosettaOutputsBase_fwd_hh
