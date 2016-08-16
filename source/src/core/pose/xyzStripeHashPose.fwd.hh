// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_core_pose_xyzStripeHashPose_fwd_hh
#define INCLUDED_core_pose_xyzStripeHashPose_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace pose {

enum PoseCoordPickMode {
	PoseCoordPickMode_NBR,
	PoseCoordPickMode_CA,
	PoseCoordPickMode_CB,
	PoseCoordPickMode_CB_else_CA,
	PoseCoordPickMode_BB,
	PoseCoordPickMode_N_CA_C,
	PoseCoordPickMode_N_CA_C_CB,
	PoseCoordPickMode_N_C_O,
	PoseCoordPickMode_BNP,
	PoseCoordPickMode_HVY,
	PoseCoordPickMode_HVY_IF_NP,
	PoseCoordPickMode_ALL,
	PoseCoordPickMode_CBorCA,
	PoseCoordPickMode_NUL
};

class xyzStripeHashPose;
typedef utility::pointer::shared_ptr< xyzStripeHashPose > xyzStripeHashPoseOP;
typedef utility::pointer::shared_ptr< xyzStripeHashPose const > xyzStripeHashPoseCOP;
typedef utility::pointer::weak_ptr< xyzStripeHashPose const > xyzStripeHashPoseCAP;


}
}

#endif
