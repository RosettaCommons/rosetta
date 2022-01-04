// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CycpepRigidBodyPermutationMover.fwd.hh
/// @brief Forward declarations for a mover that takes a cyclic peptide and alters its position (rigid-body transform), superimposing
/// it on a cyclic permutation or reverse cyclic permutation of itself.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_cyclic_peptide_CycpepRigidBodyPermutationMover_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_CycpepRigidBodyPermutationMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace cyclic_peptide {

class CycpepRigidBodyPermutationMover;

using CycpepRigidBodyPermutationMoverOP = utility::pointer::shared_ptr< CycpepRigidBodyPermutationMover >;
using CycpepRigidBodyPermutationMoverCOP = utility::pointer::shared_ptr< CycpepRigidBodyPermutationMover const >;

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_CycpepRigidBodyPermutationMover_fwd_hh
