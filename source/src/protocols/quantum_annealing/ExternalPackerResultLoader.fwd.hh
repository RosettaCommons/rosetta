// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/quantum_annealing/ExternalPackerResultLoader.fwd.hh
/// @brief A mover that reads the result of an externally-called packer, and applies it to a pose to generate a structure.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_quantum_annealing_ExternalPackerResultLoader_fwd_hh
#define INCLUDED_protocols_quantum_annealing_ExternalPackerResultLoader_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace quantum_annealing {

class ExternalPackerResultLoader;

typedef utility::pointer::shared_ptr< ExternalPackerResultLoader > ExternalPackerResultLoaderOP;
typedef utility::pointer::shared_ptr< ExternalPackerResultLoader const > ExternalPackerResultLoaderCOP;

} //protocols
} //quantum_annealing

#endif //INCLUDED_protocols_quantum_annealing_ExternalPackerResultLoader_fwd_hh
