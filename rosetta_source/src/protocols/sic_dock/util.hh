// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols
#ifndef INCLUDED_protocols_sic_dock_util_hh
#define INCLUDED_protocols_sic_dock_util_hh

#include <core/types.hh>
#include <core/kinematics/Stub.hh>
#include <protocols/sic_dock/types.hh>
#include <protocols/sic_dock/RigidScore.fwd.hh>
#include <protocols/sic_dock/SICFast.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>
//#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>


namespace protocols {
namespace sic_dock {

int
neighbor_count(
	core::pose::Pose const & pose,
	int ires,
	double distance_threshold=10.0
);

double
slide_into_contact_and_score(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	core::kinematics::Stub                & xa,
	core::kinematics::Stub          const & xb,
	numeric::xyzVector<core::Real>  const & ori,
	core::Real                            & score
);

core::pose::Pose const &                   pose_with_most_CBs( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );
bool                                       pose1_has_most_CBs( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );
core::Size                                          count_CBs( core::pose::Pose const & pose );
core::id::AtomID_Map<double>          cb_weight_map_from_pose( core::pose::Pose const & pose );
utility::vector1<numeric::xyzVector<core::Real> > get_CB_Vecs( core::pose::Pose const & pose );
utility::vector1<core::Real>             cb_weights_from_pose( core::pose::Pose const & pose );

void
xform_pose(
	core::pose::Pose & pose,
	core::kinematics::Stub const & s,
	core::Size sres=1,
	core::Size eres=0
);
void
xform_pose_rev(
	core::pose::Pose & pose,
	core::kinematics::Stub const & s,
	core::Size sres=1,
	core::Size eres=0
);

utility::vector1<core::Size> range(core::Size beg, core::Size end);
Vec3
get_leap_lower_stub(
	core::pose::Pose const & pose,
	core::Size ir
);

int flood_fill3D(int i, int j, int k, ObjexxFCL::FArray3D<double> & grid, double t);

// void termini_exposed(core::pose::Pose const & pose, bool & ntgood, bool & ctgood );

inline core::kinematics::Stub multstubs(core::kinematics::Stub const & a, core::kinematics::Stub const & b){
	return core::kinematics::Stub( a.M*b.M, a.M*b.v+a.v );
}
inline core::kinematics::Stub invstub(core::kinematics::Stub const & a){
	numeric::xyzMatrix<core::Real> const MR = a.M.transposed();
	return core::kinematics::Stub( MR, MR * -a.v );
}


} // sic_dock
} // protocols

#endif // INCLUDED_protocols_sic_dock_util_HH
