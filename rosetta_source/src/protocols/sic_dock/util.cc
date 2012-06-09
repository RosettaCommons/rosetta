// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <core/pose/util.hh>
#include <numeric/xyz.functions.hh>

namespace protocols {
namespace sic_dock {

using core::Size;
using core::Real;
using numeric::min;
using core::id::AtomID;
typedef numeric::xyzVector<core::Real> Vec;


int neighbor_count(core::pose::Pose const &pose, int ires, double distance_threshold=10.0) {
	core::conformation::Residue const resi( pose.residue( ires ) );
	Size resi_neighbors( 0 );
	for(Size jres = 1; jres <= pose.n_residue(); ++jres) {
		core::conformation::Residue const resj( pose.residue( jres ) );
		double const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
		if( distance <= distance_threshold ){
			++resi_neighbors;
		}
	}
	return resi_neighbors;
}


double
slide_into_contact_and_score(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	core::kinematics::Stub                & xa,
	core::kinematics::Stub          const & xb,
	numeric::xyzVector<core::Real>  const & ori,
	core::Real                            & score
){
	double d = sic.slide_into_contact(xa,xb,ori);
	xa.v += d*ori;
	if(score != -12345.0) score = sfxn.score( xa, xb );
	return d;
}

core::id::AtomID_Map<double>
cb_weight_map_from_pose(
	core::pose::Pose const & pose
){
	core::id::AtomID_Map<double> amap;
	core::pose::initialize_atomid_map(amap,pose,-1.0);
	for(Size i = 1; i <= pose.n_residue(); ++i){
		if(pose.residue(i).has("CB")) {
			double cbw = numeric::min(1.0,(double)neighbor_count(pose,i)/20.0);
			amap[AtomID(pose.residue(i).atom_index("CB"),i)] = cbw;
		}
	}
	return amap;
}

utility::vector1<core::Real>
cb_weights_from_pose(
	core::pose::Pose const & pose
){
	utility::vector1<core::Real> wts;
	for(Size i = 1; i <= pose.n_residue(); ++i){
		if(pose.residue(i).has("CB")) {
			double cbw = numeric::min(1.0,(double)neighbor_count(pose,i)/20.0);
			wts.push_back(cbw);
		}
	}
	return wts;
}

core::Size
count_CBs(
	core::pose::Pose const & pose
){
	core::Size cbcount;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) cbcount += pose.residue(ir).has("CB");
	return cbcount;
}

core::pose::Pose const &
pose_with_most_CBs(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2
){
	bool most_is_one = count_CBs(pose1) > count_CBs(pose2);
	return most_is_one? pose1: pose2;
}

bool
pose1_has_most_CBs(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2
){
	return count_CBs(pose1) > count_CBs(pose2);
}

utility::vector1<numeric::xyzVector<core::Real> >
get_CB_Vecs(
	core::pose::Pose const & pose
){
	utility::vector1<numeric::xyzVector<core::Real> > CBs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		if( pose.residue(ir).has("CB") ){
			CBs.push_back( pose.residue(ir).xyz("CB") );
		}
	}
	return CBs;	
}

void
xform_pose( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres=1, Size eres=0 ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
    }
  }
}

} // sic_dock
} // protocols
