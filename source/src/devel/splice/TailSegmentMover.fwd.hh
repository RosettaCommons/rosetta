// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief  Splice forward declarations header

#ifndef INCLUDED_devel_splice_FloppyTailMover_fwd_hh
#define INCLUDED_devel_splice_FloppyTailMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace splice {

//Forwards and OP typedefs
class TailSegmentMover;
typedef utility::pointer::owning_ptr< TailSegmentMover > TailSegmentMoverOP;
typedef utility::pointer::owning_ptr< TailSegmentMover const > TailSegmentMoverCOP;

} //splice
} //devel

#endif //INCLUDED_devel_splice_SpliceSegment_fwd_hh
