// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_protocols_sic_dock_loophash_util_hh
#define INCLUDED_protocols_sic_dock_loophash_util_hh

#include <core/types.hh>
#include <core/kinematics/Stub.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>
#include <protocols/sic_dock/types.hh>
#include <protocols/sic_dock/RigidScore.fwd.hh>
#include <protocols/sic_dock/SICFast.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>

//#include <numeric/xyzVector.fwd.hh>

namespace protocols {
namespace sic_dock {



Vec3
get_leap_upper_stub(
	core::pose::Pose const & pose,
	platform::Size ir
);
Xform vec3_to_stub(Vec3 const & v3);
Xform vec3_to_stub(Xform const & xform, Vec3 const & v3);

void
get_termini_from_pose(
	core::pose::Pose const & pose,
	platform::Size ir,
	TermInfo & lowers,
	TermInfo & uppers
);
void
get_termini_from_pose(
	core::pose::Pose const & pose,
	TermInfo & lowers,
	TermInfo & uppers
);
numeric::geometry::hashing::Real6
get_leap_6dof(
	Xform const & lower,
	Xform const & upper
);
platform::Size
count_linkers(
	Xform const & lower,
	Xform const & upper,
	protocols::loophash::LoopHashLibraryOP loop_hash_library,
	utility::vector1<platform::Size> const & loopsizes,
	platform::Size radius = 0
);

platform::Size
dump_loophash_linkers(
	Xform const & lower,
	Xform const & upper,
	// core::pose::Pose const & pose1,
	// core::pose::Pose const & pose2,
	protocols::loophash::LoopHashLibraryOP loop_hash_library,
	utility::vector1<platform::Size> const & loopsizes,
	platform::Size radius = 0,
	std::string const & outtag = ""
);

platform::Real
linker_count2score(
	platform::Size count
);

} // sic_dock
} // protocols

#endif // INCLUDED_protocols_sic_dock_util_HH
