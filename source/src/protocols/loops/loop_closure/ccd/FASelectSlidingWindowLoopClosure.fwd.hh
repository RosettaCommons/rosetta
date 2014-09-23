// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/ccd/FASelectSlidingWindowLoopClosure.fwd.hh
/// @brief  fwd headers for FASelectSlidingWindowLoopClosure
/// @author Brian Weitzner

#ifndef INCLUDED_protocols_loops_loop_closure_ccd_FASelectSlidingWindowLoopClosure_fwd_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_FASelectSlidingWindowLoopClosure_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

// Forward
class FASelectSlidingWindowLoopClosure;

// Types
typedef  utility::pointer::shared_ptr< FASelectSlidingWindowLoopClosure >  FASelectSlidingWindowLoopClosureOP;
typedef  utility::pointer::shared_ptr< FASelectSlidingWindowLoopClosure const >  FASelectSlidingWindowLoopClosureCOP;

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_FASelectSlidingWindowLoopClosure_fwd_hh
