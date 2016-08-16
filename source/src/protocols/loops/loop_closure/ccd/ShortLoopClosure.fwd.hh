// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_closure/ShortLoopClosure.fwd.hh
/// @brief  fwd headers for ShortLoopClosure
/// @author Brian Weitzner

#ifndef INCLUDED_protocols_loops_loop_closure_ccd_ShortLoopClosure_fwd_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_ShortLoopClosure_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

// Forward
class ShortLoopClosure;

// Types
typedef  utility::pointer::shared_ptr< ShortLoopClosure >  ShortLoopClosureOP;
typedef  utility::pointer::shared_ptr< ShortLoopClosure const >  ShortLoopClosureCOP;

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_ShortLoopClosure_fwd_hh
