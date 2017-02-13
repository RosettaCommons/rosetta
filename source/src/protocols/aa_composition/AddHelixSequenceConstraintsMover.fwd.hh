// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddHelixSequenceConstraintsMover.fwd.hh
/// @brief This mover adds sequence constraints to the ends of each helix, requiring at least one positively-charged residue in the three C-terminal residues, and at least one negatively-charged resiude in the three N-terminal residues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMover_fwd_hh
#define INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace aa_composition {

class AddHelixSequenceConstraintsMover;

typedef utility::pointer::shared_ptr< AddHelixSequenceConstraintsMover > AddHelixSequenceConstraintsMoverOP;
typedef utility::pointer::shared_ptr< AddHelixSequenceConstraintsMover const > AddHelixSequenceConstraintsMoverCOP;

} //protocols
} //aa_composition

#endif //INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMover_fwd_hh
