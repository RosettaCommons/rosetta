// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/sic_dock/util.hh>
#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <core/pose/util.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray3D.hh>

namespace protocols {
namespace sic_dock {

using platform::Size;
using platform::Real;
using numeric::min;
using core::id::AtomID;
using std::cout;
using std::endl;
typedef platform::Real Real;
typedef platform::Size Size;
typedef core::pose::Pose Pose;
typedef core::kinematics::Stub Stub;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;
typedef utility::vector1<Real> Reals;
typedef utility::vector1<Size> Sizes;
typedef utility::vector1<Stub> Stubs;
typedef utility::vector1<RigidScoreCOP> Scores;


int neighbor_count(core::pose::Pose const &pose, int ires, double distance_threshold) {
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

Real
cb_weight(core::pose::Pose const &pose, Size ires, Real distance_threshold) {
	Real wt = numeric::min(1.0,(double)neighbor_count(pose,ires,distance_threshold)/20.0);
	if(pose.secstruct(ires)=='L') wt = wt / 3.0; //TODO make option somehow
	return wt;
}

double
slide_into_contact_and_score(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	numeric::xyzTransform<core::Real>         & xa,
	numeric::xyzTransform<core::Real>   const & xb,
	numeric::xyzVector<platform::Real>  const & ori,
	platform::Real                            & score
){
	Stub sa(xa.R,xa.t), sb(xb.R,xb.t);
	double d = sic.slide_into_contact(sa,sb,ori);
	sa.v += d*ori;
	xa.t += d*ori;
	if(score != -12345.0) score = sfxn.score( sa, sb );
	return d;	
}

double
slide_into_contact_and_score(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	core::kinematics::Stub                & xa,
	core::kinematics::Stub          const & xb,
	numeric::xyzVector<platform::Real>  const & ori,
	platform::Real                            & score
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
			amap[AtomID(pose.residue(i).atom_index("CB"),i)] = cb_weight(pose,i);
		}
	}
	return amap;
}

utility::vector1<platform::Real>
cb_weights_from_pose(
	core::pose::Pose const & pose
){
	utility::vector1<platform::Real> wts;
	for(Size i = 1; i <= pose.n_residue(); ++i){
		if(pose.residue(i).has("CB")) {
			wts.push_back( cb_weight(pose,i) );
		}
	}
	return wts;
}

platform::Size
count_CBs(
	core::pose::Pose const & pose
){
	platform::Size cbcount;
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

utility::vector1<numeric::xyzVector<platform::Real> >
get_CB_Vecs(
	core::pose::Pose const & pose
){
	utility::vector1<numeric::xyzVector<platform::Real> > CBs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		if( pose.residue(ir).has("CB") ){
			CBs.push_back( pose.residue(ir).xyz("CB") );
		}
	}
	return CBs;	
}

void
xform_pose( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres, Size eres ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
    }
  }
}

void
xform_pose_rev( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres, Size eres ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s.global2local(pose.xyz(aid)) );
    }
  }
}

void xform_pose( core::pose::Pose & pose, numeric::xyzTransform<core::Real> const & s, Size sres, Size eres ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
}
void xform_pose_rev( core::pose::Pose & pose, numeric::xyzTransform<core::Real> const & s ) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, ~s * pose.xyz(aid) );
    }
  }
}



utility::vector1<platform::Size> range(platform::Size beg, platform::Size end){
	utility::vector1<platform::Size> v;
	for(platform::Size i = beg; i < end; ++i) v.push_back(i);
	return v;
}

int flood_fill3D(int i, int j, int k, ObjexxFCL::FArray3D<double> & grid, double t) {
	if( grid(i,j,k) <= t ) return 0;
	grid(i,j,k) = t;
	int nmark = 1;
	if(i>1                ) nmark += flood_fill3D(i-1,j  ,k  ,grid,t);
	if(i<(int)grid.size1()) nmark += flood_fill3D(i+1,j  ,k  ,grid,t);
	if(j>1                ) nmark += flood_fill3D(i  ,j-1,k  ,grid,t);
	if(j<(int)grid.size2()) nmark += flood_fill3D(i  ,j+1,k  ,grid,t);
	if(k>1                ) nmark += flood_fill3D(i  ,j  ,k-1,grid,t);
	if(k<(int)grid.size3()) nmark += flood_fill3D(i  ,j  ,k+1,grid,t);
	return nmark;
}



// void
// termini_exposed(
// 	core::pose::Pose const & pose,
// 	bool & ntgood,
// 	bool & ctgood
// ){
// 	using basic::options::option;
// 	using namespace basic::options::OptionKeys;
// 	core::id::AtomID_Map<Real> atom_sasa;
// 	core::id::AtomID_Map<bool> atom_subset;
// 	utility::vector1<Real> rsd_sasa;
// 	core::pose::initialize_atomid_map(atom_subset, pose, false);
// 	for(int i = 2; i <= (int)pose.n_residue()-1; ++i) {
// 		for(int ia = 1; ia <= (int)pose.residue(i).nheavyatoms(); ++ia) {
// 			if(pose.residue(i).atom_is_backbone(ia))
// 				atom_subset[core::id::AtomID(ia,i)] = true;
// 		}
// 	}
// 	atom_subset[core::id::AtomID(1,1)] = true;
// 	atom_subset[core::id::AtomID(3,pose.n_residue())] = true;
// 	core::scoring::calc_per_atom_sasa( pose, atom_sasa,rsd_sasa, 4.0, false, atom_subset );
// 	Real nexpose = atom_sasa[core::id::AtomID(1,        1       )] / 12.56637 / 5.44 / 5.44;
// 	Real cexpose = atom_sasa[core::id::AtomID(3,pose.n_residue())] / 12.56637 / 5.44 / 5.44;

// 	Vec nt = pose.residue(        1       ).xyz("N");
// 	Vec ct = pose.residue(pose.n_residue()).xyz("C");
// 	Real nang = angle_degrees(nt,Vec(0,0,0),Vec(nt.x(),nt.y(),0));
// 	Real cang = angle_degrees(ct,Vec(0,0,0),Vec(ct.x(),ct.y(),0));
// 	ntgood = nexpose > option[sicdock::term_min_expose]() && nang < option[sicdock::term_max_angle]();
// 	ctgood = cexpose > option[sicdock::term_min_expose]() && cang < option[sicdock::term_max_angle]();
// 	// platform::Real nnt=0.0,nct=0.0,gnt=0.0,gct=0.0;
// 	// for(int ir=1; ir<=pose.n_residue(); ++ir) {
// 	// 	for(int ia=1; ia<=5; ++ia) {
// 	// 		Vec x = pose.residue(ir).xyz(ia);
// 	// 		if(angle_degrees(x,Vec(0,0,0),nt) < 15.0 &&  ) {
// 	// 			nnt += 1.0;
// 	// 			if( nt.normalized().dot(x) < nt.length() )
// 	// 				gnt += 1.0;
// 	// 		}
// 	// 		if(angle_degrees(x,Vec(0,0,0),ct) < 15.0 ) {
// 	// 			nct += 1.0;
// 	// 			if( ct.normalized().dot(x) < ct.length() )
// 	// 				gct += 1.0;
// 	// 		}
// 	// 	}
// 	// }
// }



} // sic_dock
} // protocols
