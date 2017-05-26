// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_setup_FullModelInfoSetupFromCommandLine_HH
#define INCLUDED_protocols_stepwise_setup_FullModelInfoSetupFromCommandLine_HH

#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>
#include <map>

namespace protocols {
namespace stepwise {
namespace setup {

void
initialize_native_and_align_pose( core::pose::PoseOP & native_pose,
	core::pose::PoseOP & align_pose,
	core::chemical::ResidueTypeSetCAP rsd_set,
	core::pose::PoseCOP start_pose );


} //setup
} //stepwise
} //protocols

#endif
