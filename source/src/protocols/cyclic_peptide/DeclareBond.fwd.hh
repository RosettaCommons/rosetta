// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide/DeclareBond.fwd.hh
/// @brief  Defines owning pointers for DeclareBond mover class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_DeclareBond_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_DeclareBond_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace cyclic_peptide {

class DeclareBond; // fwd declaration
typedef utility::pointer::shared_ptr< DeclareBond > DeclareBondOP;
typedef utility::pointer::shared_ptr< DeclareBond const > DeclareBondCOP;
typedef utility::vector1<DeclareBondOP> DeclareBondOPs;
typedef utility::vector1<DeclareBondCOP> DeclareBondCOPs;

} // namespace cyclic_peptide
} // namespace protocols

#endif // INCLUDED_protocols_cyclic_peptide_DeclareBond_fwd_hh
