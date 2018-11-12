// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddMHCEpitopeConstraintMover.fwd.hh
/// @brief Forward declarations for the AddMHCEpitopeConstraintMover.
/// @details Assigns an MHCEpitopeConstraint to a pose, initializing it from a file.
/// @author Chris Bailey-Kellogg (cbk@cs.dartmouth.edu), based on Vikram Mulligan's MHCEpitopeConstraint

#ifndef INCLUDED_protocols_aa_composition_AddMHCEpitopeConstraintMover_fwd_hh
#define INCLUDED_protocols_aa_composition_AddMHCEpitopeConstraintMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace aa_composition {

class AddMHCEpitopeConstraintMover;
typedef utility::pointer::shared_ptr< AddMHCEpitopeConstraintMover > AddMHCEpitopeConstraintMoverOP;
typedef utility::pointer::shared_ptr< AddMHCEpitopeConstraintMover const > AddMHCEpitopeConstraintMoverCOP;

} // aa_composition
} // protocols

#endif
