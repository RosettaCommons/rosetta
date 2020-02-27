// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.fwd.hh
/// @brief Loop closure mover that uses linearization of the kinematics to adjust all unlocked backbone torsion angles of a loop at the same time
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

#ifndef INCLUDED_protocols_loops_loop_closure_jacobi_JacobiLoopClosureMover_fwd_hh
#define INCLUDED_protocols_loops_loop_closure_jacobi_JacobiLoopClosureMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace loops {
namespace loop_closure {
namespace jacobi {

class JacobiLoopClosureMover;

using JacobiLoopClosureMoverOP = utility::pointer::shared_ptr< JacobiLoopClosureMover >;
using JacobiLoopClosureMoverCOP = utility::pointer::shared_ptr< JacobiLoopClosureMover const >;

} //jacobi
} //loop_closure
} //loops
} //protocols

#endif //INCLUDED_protocols_loops_loop_closure_jacobi_JacobiLoopClosureMover_fwd_hh
