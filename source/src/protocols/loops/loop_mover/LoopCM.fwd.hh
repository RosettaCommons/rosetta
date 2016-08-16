// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_mover/LoopCM.fwd.hh
/// @brief  LoopCM forward declarations
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_loops_loop_mover_LoopCM_fwd_hh
#define INCLUDED_protocols_loops_loop_mover_LoopCM_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_mover {

// Forward
class LoopCM;

typedef utility::pointer::shared_ptr< LoopCM > LoopCMOP;
typedef utility::pointer::shared_ptr< LoopCM const > LoopCMCOP;

typedef utility::pointer::weak_ptr< LoopCM > LoopCMAP;
typedef utility::pointer::weak_ptr< LoopCM const > LoopCMCAP;

} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_LoopMover_FWD_HH

