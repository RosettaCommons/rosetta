// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#ifndef INCLUDED_protocols_sic_dock_util_hh
#define INCLUDED_protocols_sic_dock_util_hh

#include <core/types.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/sic_dock/types.hh>
#include <protocols/sic_dock/RigidScore.fwd.hh>
#include <protocols/sic_dock/SICFast.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <numeric/xyzTransform.hh>
#include <core/kinematics/Stub.fwd.hh>

#include <numeric/xyzVector.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/xyzTriple.hh>
#include <numeric/Quaternion.hh>


namespace protocols {
namespace sic_dock {

core::Real get_rg(core::pose::Pose const & p);

int
neighbor_count(
	core::pose::Pose const & pose,
	int ires,
	double distance_threshold=10.0
);

core::Real
cb_weight(
	core::pose::Pose const &pose,
	core::Size ires,
	core::Real distance_threshold=10.0
);

void make_Cx(core::pose::Pose & pose, int N, numeric::xyzVector<core::Real> axis=numeric::xyzVector<core::Real>(0,0,1) );

double
slide_into_contact_and_score(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	Xform                                 & xa,
	Xform                           const & xb,
	numeric::xyzVector<core::Real>  const & ori,
	core::Real                            & score
);

double
slide_into_contact_and_score_DEPRICATED(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	core::kinematics::Stub                & xa,
	core::kinematics::Stub          const & xb,
	numeric::xyzVector<core::Real> const  & ori,
	core::Real                            & score
);

core::pose::Pose const &                   pose_with_most_CBs( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );
bool                                       pose1_has_most_CBs( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

core::Size                                          count_CBs( core::pose::Pose const & pose );

core::id::AtomID_Map<double>          cb_weight_map_from_pose( core::pose::Pose const & pose );

utility::vector1<numeric::xyzVector<core::Real> > get_CB_Vecs_from_pose( core::pose::Pose const & pose );
utility::vector1<numeric::xyzVector<core::Real> > get_CB_Vecs_from_map ( core::pose::Pose const & pose, core::id::AtomID_Map<core::Real> const & map );

utility::vector1<                   core::Real>    cb_weights_from_pose( core::pose::Pose const & pose );
utility::vector1<                   core::Real>    cb_weights_from_map ( core::pose::Pose const & pose, core::id::AtomID_Map<core::Real> const & map );

utility::vector1<core::Size> range(core::Size beg, core::Size end);
Vec3
get_leap_lower_stub(
	core::pose::Pose const & pose,
	core::Size ir
);

int flood_fill3D(int i, int j, int k, ObjexxFCL::FArray3D<double> & grid, double t);

// void termini_exposed(core::pose::Pose const & pose, bool & ntgood, bool & ctgood );

inline Xform multstubs(Xform const & a, Xform const & b){
	return Xform( a.R*b.R, a.R*b.t+a.t );
}
inline Xform invstub(Xform const & a){
	numeric::xyzMatrix<core::Real> const MR = a.R.transposed();
	return Xform( MR, MR * -a.t );
}

bool residue_is_floppy(core::pose::Pose const & pose, core::Size const ir, core::Real const ttrim_cut=1.0, core::Size const nfold=1);

void auto_trim_floppy_termini(core::pose::Pose & pose, core::Real const ttrim_cut=1.0, core::Size const nfold=1);

numeric::xyzVector<core::Real> center_of_geom(core::pose::Pose const & pose, core::Size str=1, core::Size end=0);
void dump_points_pdb(utility::vector1<numeric::xyzVector<core::Real> > const & p, std::string fn);
void dump_points_pdb(utility::vector1<numeric::xyzVector<core::Real> > const & p, numeric::xyzVector<core::Real>  t, std::string fn);
void trans_pose(  core::pose::Pose & pose, numeric::xyzVector<core::Real>  const & trans, core::Size start=1, core::Size end=0 );
void rot_pose  (  core::pose::Pose & pose, numeric::xyzMatrix<core::Real> const & rot, core::Size start=1, core::Size end=0 );
void rot_pose  (  core::pose::Pose & pose, numeric::xyzMatrix<core::Real> const & rot, numeric::xyzVector<core::Real>  const & cen, core::Size start=1, core::Size end=0 );
void rot_pose  (  core::pose::Pose & pose, numeric::xyzVector<core::Real> const & axis, core::Real const & ang, core::Size start=1, core::Size end=0 );
void rot_pose  (  core::pose::Pose & pose, numeric::xyzVector<core::Real> const & axis, core::Real const & ang, numeric::xyzVector<core::Real>  const & cen, core::Size start=1, core::Size end=0 );
void alignaxis ( core::pose::Pose & pose, numeric::xyzVector<core::Real>  newaxis, numeric::xyzVector<core::Real>  oldaxis, numeric::xyzVector<core::Real>  cen = numeric::xyzVector<core::Real> (0,0,0) );
numeric::xyzTransform<core::Real> alignaxis_xform (numeric::xyzVector<core::Real>  newaxis, numeric::xyzVector<core::Real>  oldaxis, numeric::xyzVector<core::Real>  cen = numeric::xyzVector<core::Real> (0,0,0) );


inline numeric::xyzVector<core::Real>  projperp(numeric::xyzVector<core::Real>  const & u, numeric::xyzVector<core::Real>  const & v) {  return v - projection_matrix(u)*v; }

void xform_pose ( core::pose::Pose & pose,            core::kinematics::Stub const & s, core::Size sres=1, core::Size eres=0 );
void xform_pose ( core::pose::Pose & pose, numeric::xyzTransform<core::Real> const & s, core::Size sres=1, core::Size eres=0 );
void xform_pose_rev (core::pose::Pose & pose, core::kinematics::Stub const & s);
void xform_pose_rev (core::pose::Pose & pose, numeric::xyzTransform<core::Real> const & s);

core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::Residue const & fixd_resi);

template<typename T> inline T sqr(T x) { return x*x; }

std::string KMGT(double const & x, int const & w, int const & d);

} // sic_dock
} // protocols

#endif // INCLUDED_protocols_sic_dock_util_HH
