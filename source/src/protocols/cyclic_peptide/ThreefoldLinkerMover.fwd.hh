// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/ThreefoldLinkerMover.fwd.hh
/// @brief This mover links three cysteine residues with a three-way cross-linker.  It adds the crosslinker,
/// sets up constraints, optionally packs and energy-mimizes it into place (packing/minimizing only the crosslinker and
/// the side-chains to which it connects), andthen optionally relaxes the whole structure.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_protocols_cyclic_peptide_ThreefoldLinkerMover_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_ThreefoldLinkerMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace cyclic_peptide {

class ThreefoldLinkerMover;

typedef utility::pointer::shared_ptr< ThreefoldLinkerMover > ThreefoldLinkerMoverOP;
typedef utility::pointer::shared_ptr< ThreefoldLinkerMover const > ThreefoldLinkerMoverCOP;

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_ThreefoldLinkerMover_fwd_hh
