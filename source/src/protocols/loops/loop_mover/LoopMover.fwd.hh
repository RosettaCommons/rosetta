// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_mover/LoopMover.fwd.hh
/// @brief  LoopMover base classes
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loops_loop_mover_LoopMover_fwd_hh
#define INCLUDED_protocols_loops_loop_mover_LoopMover_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_mover {

// Forward
class LoopMover;

typedef utility::pointer::shared_ptr< LoopMover > LoopMoverOP;
typedef utility::pointer::shared_ptr< LoopMover const > LoopMoverCOP;

typedef utility::pointer::weak_ptr< LoopMover > LoopMoverAP;
typedef utility::pointer::weak_ptr< LoopMover const > LoopMoverCAP;

class MultiLoopMover;

typedef utility::pointer::shared_ptr< MultiLoopMover > MultiLoopMoverOP;
typedef utility::pointer::shared_ptr< MultiLoopMover const > MultiLoopMoverCOP;

class SingleLoopMover;

typedef utility::pointer::shared_ptr< SingleLoopMover > SingleLoopMoverOP;
typedef utility::pointer::shared_ptr< SingleLoopMover const > SingleLoopMoverCOP;

} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_LoopMover_FWD_HH

