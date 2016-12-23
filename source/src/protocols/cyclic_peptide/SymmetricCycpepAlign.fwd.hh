// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/SymmetricCycpepAlign.fwd.hh
/// @brief Given a quasi-symmetric cyclic peptide, this mover aligns the peptide so that the cyclic symmetry axis lies along the Z-axis and the centre of mass is at the origin.
/// It then optionally removes all but one symmetry repeat, so that true symmetry may be set up with the SetupForSymmetry mover.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_SymmetricCycpepAlign_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_SymmetricCycpepAlign_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace cyclic_peptide {

class SymmetricCycpepAlign;

typedef utility::pointer::shared_ptr< SymmetricCycpepAlign > SymmetricCycpepAlignOP;
typedef utility::pointer::shared_ptr< SymmetricCycpepAlign const > SymmetricCycpepAlignCOP;

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_SymmetricCycpepAlign_fwd_hh
