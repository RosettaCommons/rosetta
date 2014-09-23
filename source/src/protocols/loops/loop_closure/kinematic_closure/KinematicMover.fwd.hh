// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh
/// @brief  KinematicMover forward declarations header
/// @author Steven Lewis (smlewi@unc.edu)


#ifndef INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicMover_fwd_hh
#define INCLUDED_protocols_loops_loop_closure_kinematic_closure_KinematicMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

//Forwards and OP typedefs
class KinematicMover;
typedef utility::pointer::shared_ptr< KinematicMover > KinematicMoverOP;
typedef utility::pointer::shared_ptr< KinematicMover const > KinematicMoverCOP;

typedef utility::pointer::weak_ptr< KinematicMover const > KinematicMoverCAP;

} // namespace kinematic_closure
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_KinematicMover_FWD_HH
