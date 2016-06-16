// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cyclic_peptide/FlipChiralityMover.fwd.hh
/// @brief  Defines owning pointers for FlipChiralityMover mover class.
/// @author Parisa Hosseinzadeh (parisah@uw.edu) and Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_FlipChiralityMover_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_FlipChiralityMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace cyclic_peptide {

class FlipChiralityMover; // fwd declaration
typedef utility::pointer::shared_ptr< FlipChiralityMover > FlipChiralityMoverOP;
typedef utility::pointer::shared_ptr< FlipChiralityMover const > FlipChiralityMoverCOP;
typedef utility::vector1<FlipChiralityMoverOP> FlipChiralityMoverOPs;
typedef utility::vector1<FlipChiralityMoverCOP> FlipChiralityMoverCOPs;

} // namespace cyclic_peptide
} // namespace protocols

#endif // INCLUDED_protocols_cyclic_peptide_FlipChiralityMover_fwd_hh
