// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_xyzStripeHashPose_hh
#define INCLUDED_protocols_sic_dock_xyzStripeHashPose_hh

#include <protocols/sic_dock/xyzStripeHashPose.fwd.hh>
#include <protocols/sic_dock/types.hh>
#include <numeric/geometry/hashing/xyzStripeHash.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>

#include <platform/types.hh>

namespace protocols {
namespace sic_dock {


class xyzStripeHashPose : public numeric::geometry::hashing::xyzStripeHash {
public:
	xyzStripeHashPose(
		double radius=0.0
	);
	xyzStripeHashPose(
		core::pose::Pose const & p,
		PoseCoordPickMode m = BB,
		double radius = 0.0
	);
	xyzStripeHashPose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<double> const & amap,
		double radius = 0.0
	);
	void
	add_pose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<double> const amap
	);
	void
	add_pose(
		core::pose::Pose const & p,
		PoseCoordPickMode m = BB
	);
	void
	init_posehash();
private:
	utility::vector1<numeric::geometry::hashing::Ball> balls_;
	bool initialized_;
};


} // namespace sic_dock
} // namespace protocols

#endif
