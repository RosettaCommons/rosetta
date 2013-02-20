// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/loophash/LoopHashLoopClosureMover.fwd.hh
/// @brief  LoopHashLoopClosureMover class forward delcaration
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_protocols_loops_loop_closure_loophash_LoopHashLoopClosureMover_HH
#define INCLUDED_protocols_loops_loop_closure_loophash_LoopHashLoopClosureMover_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace loophash {
class  LoopHashLoopClosureMover;
typedef utility::pointer::owning_ptr< LoopHashLoopClosureMover > LoopHashLoopClosureMoverOP;
typedef utility::pointer::owning_ptr< LoopHashLoopClosureMover const > LoopHashLoopClosureMoverCOP;
typedef utility::pointer::access_ptr< LoopHashLoopClosureMover > LoopHashLoopClosureMoverAP;
typedef utility::pointer::access_ptr< LoopHashLoopClosureMover const > LoopHashLoopClosureMoverCAP;

}
}
}
}
#endif
