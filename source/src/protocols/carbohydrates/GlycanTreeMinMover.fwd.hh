// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeMinMover.fwd.hh
/// @brief A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_carbohydrates_GlycanTreeMinMover_fwd_hh
#define INCLUDED_protocols_carbohydrates_GlycanTreeMinMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace carbohydrates {

class GlycanTreeMinMover;

typedef utility::pointer::shared_ptr< GlycanTreeMinMover > GlycanTreeMinMoverOP;
typedef utility::pointer::shared_ptr< GlycanTreeMinMover const > GlycanTreeMinMoverCOP;



} //protocols
} //carbohydrates


#endif //INCLUDED_protocols_carbohydrates_GlycanTreeMinMover_fwd_hh





