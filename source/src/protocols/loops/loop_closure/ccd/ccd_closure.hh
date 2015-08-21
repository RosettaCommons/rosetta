// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/ccd.ccd_closure.hh
/// @brief  Method declarations for cyclic coordinate descent loop closure.
/// @author Phil Bradley
/// @author Brian Weitzner
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_ccd_closure_HH
#define INCLUDED_protocols_loops_loop_closure_ccd_ccd_closure_HH


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

std::pair< core::Real, core::Real > get_deviation(
	core::pose::Pose const & pose, core::uint const cutpoint );

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_ccd_closure_HH
