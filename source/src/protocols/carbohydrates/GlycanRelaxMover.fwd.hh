// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanRelaxMover.fwd.hh
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_carbohydrates_GlycanRelaxMover_fwd_hh
#define INCLUDED_protocols_carbohydrates_GlycanRelaxMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace carbohydrates {

class GlycanRelaxMover;

typedef utility::pointer::shared_ptr< GlycanRelaxMover > GlycanRelaxMoverOP;
typedef utility::pointer::shared_ptr< GlycanRelaxMover const > GlycanRelaxMoverCOP;



} //protocols
} //carbohydrates


#endif //INCLUDED_protocols/carbohydrates_GlycanRelaxMover_fwd_hh





