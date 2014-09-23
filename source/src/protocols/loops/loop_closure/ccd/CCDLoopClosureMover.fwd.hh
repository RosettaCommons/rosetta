// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/loops/loop_closure/ccd/CCDLoopClosureMover.fwd.hh
/// @brief   Foward declarations for CCDLoopClosureMover
/// @author  Phil Bradley
/// @author  Oliver Lange
/// @author  Brian Weitzner
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    This file is the result of a refactor of code written by Phil and later wrapped in a Mover by Oliver.

#ifndef INCLUDED_protocols_loops_loop_closure_ccd_CCDLoopClosureMover_FWD_HH
#define INCLUDED_protocols_loops_loop_closure_ccd_CCDLoopClosureMover_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

/// @brief A Mover that performs a cyclic coordination descent loop closure on a Pose.
class CCDLoopClosureMover;

// Types
typedef  utility::pointer::shared_ptr< CCDLoopClosureMover >  CCDLoopClosureMoverOP;
typedef  utility::pointer::shared_ptr< CCDLoopClosureMover const >  CCDLoopClosureMoverCOP;

// I'm not sure what the heck this is yet; we'll just leave it here for now.... ~Labonte
class CcdMover;
typedef utility::pointer::shared_ptr< CcdMover > CcdMoverOP;
typedef utility::pointer::shared_ptr< CcdMover const > CcdMoverCOP;

}  // namespace ccd
}  // namespace loop_closure
}  // namespace loops
}  // namespace protocols

#endif  // INCLUDED_protocols_loops_loop_closure_ccd_CCDLoopClosureMover_FWD_HH
