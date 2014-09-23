// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_mover/refine/LoopMover_Backrub.fwd.hh
/// @brief  LoopMover_Backrub forward declaration
/// @author Brian Weitzner


#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCD_fwd_hh
#define INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCD_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

// Forward
class LoopMover_Refine_CCD;

typedef utility::pointer::shared_ptr< LoopMover_Refine_CCD > LoopMover_Refine_CCDOP;
typedef utility::pointer::shared_ptr< LoopMover_Refine_CCD const > LoopMover_Refine_CCDCOP;

typedef utility::pointer::weak_ptr< LoopMover_Refine_CCD > LoopMover_Refine_CCDAP;
typedef utility::pointer::weak_ptr< LoopMover_Refine_CCD const > LoopMover_Refine_CCDCAP;
} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCD_fwd_hh

