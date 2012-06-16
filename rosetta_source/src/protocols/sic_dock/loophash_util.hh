// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols
#ifndef INCLUDED_protocols_sic_dock_loophash_util_hh
#define INCLUDED_protocols_sic_dock_loophash_util_hh

#include <core/types.hh>
#include <core/kinematics/Stub.hh>
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
	core::Size ir
);
core::kinematics::Stub vec3_to_stub(Vec3 const & v3);
core::kinematics::Stub vec3_to_stub(core::kinematics::Stub const & xform, Vec3 const & v3);

void
get_termini_from_pose(
	core::pose::Pose const & pose,
	core::Size ir,
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
	core::kinematics::Stub const & lower,
	core::kinematics::Stub const & upper
);
core::Size
count_linkers(
	core::kinematics::Stub const & lower,
	core::kinematics::Stub const & upper,
	protocols::loophash::LoopHashLibraryOP loop_hash_library,
	utility::vector1<core::Size> const & loopsizes,
	core::Size radius = 0
);

core::Size
dump_loophash_linkers(
	core::kinematics::Stub const & lower,
	core::kinematics::Stub const & upper,
	// core::pose::Pose const & pose1,
	// core::pose::Pose const & pose2,
	protocols::loophash::LoopHashLibraryOP loop_hash_library,
	utility::vector1<core::Size> const & loopsizes,
	core::Size radius = 0
);

/*inline
core::kinematics::Stub
multstubs(
	core::kinematics::Stub const & a,
	core::kinematics::Stub const & b
);
inline core::kinematics::Stub invstub(core::kinematics::Stub const & a);
*/

core::Real
linker_count2score(
	core::Size count
);

} // sic_dock
} // protocols

#endif // INCLUDED_protocols_sic_dock_util_HH
