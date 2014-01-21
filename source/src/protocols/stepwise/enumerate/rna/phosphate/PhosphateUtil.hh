// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/phosphate/PhosphateUtil.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_enumerate_rna_phosphate_PhosphateUtil_HH
#define INCLUDED_protocols_stepwise_enumerate_rna_phosphate_PhosphateUtil_HH

#include <protocols/stepwise/enumerate/rna/phosphate/PhosphateMove.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace phosphate {

	void
	remove_terminal_phosphates( pose::Pose & pose );

	void
	remove_terminal_phosphates( pose::Pose & pose, utility::vector1< Size > const & res_list  );

	void
	correctly_position_five_prime_phosphate( pose::Pose & pose, Size const res );

	void
	copy_over_phosphate_variants( pose::Pose & pose,
																pose::Pose const & reference_pose,
																utility::vector1< PhosphateMove > const & phosphate_move_list );
	void
	copy_over_phosphate_variants( pose::Pose & pose,
																pose::Pose const & reference_pose,
																utility::vector1< Size > const & res_list );

} //phosphate
} //rna
} //enumerate
} //stepwise
} //protocols

#endif
