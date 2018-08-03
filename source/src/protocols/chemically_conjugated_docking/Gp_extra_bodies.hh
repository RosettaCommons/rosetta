// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/chemically_conjugated_docking/Gp_extra_bodies.hh
/// @brief contains helper quantification metrics for the original publication of the UBQ_Gp code
/// @author Steven Lewis and Hope Anderson

#ifndef INCLUDED_protocols_chemically_conjugated_docking_Gp_extra_bodies_HH
#define INCLUDED_protocols_chemically_conjugated_docking_Gp_extra_bodies_HH

// Unit Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <set>

namespace protocols {
namespace chemically_conjugated_docking {


////////////////////////////extra bodies/////////////////////////////////////////////////
/// @brief The purpose of this code is to allow for static extra things
///in the system; for its original incarnations, to allow for a RING
///domain to occlude its site on E2 and a GAP (+MG +GDP) to occlude ras
utility::vector1< core::Size > add_extra_bodies( core::pose::Pose & pose );

void pack_extra_bodies(
	utility::vector1< core::Size > const & extra_bodies_chains,
	core::pose::Pose const & pose,
	std::set< core::Size > & region);

}//chemically_conjugated_docking
}//protocols
#endif //INCLUDED_protocols_chemically_conjugated_docking_Gp_extra_bodies_HH
