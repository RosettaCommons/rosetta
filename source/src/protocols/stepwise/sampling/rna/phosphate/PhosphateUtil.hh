// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/phosphate/PhosphateUtil.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_rna_phosphate_PhosphateUtil_HH
#define INCLUDED_protocols_stepwise_sampling_rna_phosphate_PhosphateUtil_HH

#include <protocols/stepwise/sampling/rna/phosphate/PhosphateMove.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

/*
using namespace core;

Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using
*/

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace phosphate {

	void
	remove_terminal_phosphates( core::pose::Pose & pose );

	void
	remove_terminal_phosphates( core::pose::Pose & pose, utility::vector1< Size > const & res_list  );

	void
	correctly_position_five_prime_phosphate( core::pose::Pose & pose, core::Size const res );

	void
	copy_over_phosphate_variants( core::pose::Pose & pose,
								  core::pose::Pose const & reference_pose,
								  utility::vector1< PhosphateMove > const & phosphate_move_list );
	void
	copy_over_phosphate_variants( core::pose::Pose & pose,
								  core::pose::Pose const & reference_pose,
								  utility::vector1< Size > const & res_list );

} //phosphate
} //rna
} //sampling
} //stepwise
} //protocols

#endif
