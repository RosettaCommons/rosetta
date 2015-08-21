// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_pose_xyzStripeHashPose_hh
#define INCLUDED_core_pose_xyzStripeHashPose_hh

#include <core/pose/xyzStripeHashPose.fwd.hh>
#include <numeric/geometry/hashing/xyzStripeHash.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>

#include <platform/types.hh>

namespace core {
namespace pose {

core::id::AtomID_Map<platform::Real>
make_atom_map(
	core::pose::Pose const & p,
	PoseCoordPickMode m
);

class xyzStripeHashPose : public numeric::geometry::hashing::xyzStripeHash {
public:
	xyzStripeHashPose(
		platform::Real radius=0.0
	);
	xyzStripeHashPose(
		core::pose::Pose const & p,
		PoseCoordPickMode m = PoseCoordPickMode_BB,
		platform::Real radius = 0.0
	);
	xyzStripeHashPose(
		core::pose::Pose const & p,
		utility::vector1<int> const & resmap,
		PoseCoordPickMode m = PoseCoordPickMode_BB,
		platform::Real radius = 0.0
	);
	xyzStripeHashPose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<platform::Real> const & amap,
		platform::Real radius = 0.0
	);
	void
	add_pose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<platform::Real> const & amap
	);
	void
	add_pose(
		core::pose::Pose const & p,
		PoseCoordPickMode m = PoseCoordPickMode_BB
	);

	void
	init_posehash();

	static void
	extract_pose_balls(
		core::pose::Pose const & p,
		utility::vector1<numeric::geometry::hashing::Ball> & balls,
		PoseCoordPickMode m = PoseCoordPickMode_BB
	);

	static void
	extract_pose_balls(
		core::pose::Pose const & p,
		utility::vector1<numeric::geometry::hashing::Ball> & balls,
		core::id::AtomID_Map<platform::Real> const & amap
	);

private:
	utility::vector1<numeric::geometry::hashing::Ball> balls_;
	bool initialized_;
};


} // namespace pose
} // namespace core

#endif
