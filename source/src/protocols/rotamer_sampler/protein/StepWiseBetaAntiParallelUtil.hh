// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_BetaAntiparallelUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseBetaAntiparallelUtil_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseBetaAntiparallelUtil_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

//Auto Headers
#include <core/kinematics/Stub.fwd.hh>
#include <utility/vector1.hh>


using namespace::core;

namespace protocols {
namespace rotamer_sampler {
namespace protein {

	void
	do_set_xyz( pose::Pose const & pose, Size const i, pose::Pose & scratch_pose, Size const i_scratch, core::kinematics::Stub const & stub );

	void
	generate_beta_database_test();

} //protein
} //rotamer_sampler
} //protocols

#endif
