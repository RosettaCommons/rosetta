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
#include <core/kinematics/Stub.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace sic_dock {

using core::Size;
using core::Real;
using numeric::min;
using core::id::AtomID;
using std::cout;
using std::endl;
using utility::vector1;
typedef core::Real Real;
typedef core::Size Size;
typedef core::pose::Pose Pose;
typedef Xform Xform;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef vector1<Vec> Vecs;
typedef vector1<Real> Reals;
typedef vector1<Size> Sizes;
typedef numeric::Xforms Xforms;
typedef vector1<RigidScoreCOP> Scores;

static basic::Tracer TR( "protocols.sic_dock.util" );

Real get_rg(core::pose::Pose const & pose){
	Vec center_of_mass = center_of_geom(pose);
	Real rg_score = 0;
	Size nres_counted=0;
	for ( Size i = 1; i <= pose.n_residue(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		++nres_counted;
		Vec const v( pose.residue(i).nbr_atom_xyz() );
		rg_score += v.distance_squared( center_of_mass );
	}
	rg_score /= (nres_counted - 1);
	return sqrt( rg_score );
}

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

void make_Cx(core::pose::Pose & pose, int N, numeric::xyzVector<core::Real> axis ){
		core::pose::Pose tmp(pose);
		for(int inf = 2; inf <= N; ++inf){
			rot_pose( tmp, axis, 360.0/(Real)N );
			for(Size i = 1; i <= tmp.n_residue(); ++i){
				if(i==1||pose.residue(i).is_lower_terminus()||pose.residue(i).is_ligand()){
					pose.append_residue_by_jump(tmp.residue(i),1,"","",true);
				} else {
					pose.append_residue_by_bond(tmp.residue(i));
				}
			}
		}
	 }


double
slide_into_contact_and_score(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	Xform                & xa,
	Xform          const & xb,
	numeric::xyzVector<core::Real>  const & ori,
	core::Real                            & score
){
	double d = sic.slide_into_contact(xa,xb,ori);
	xa.t += d*ori;
	if(score != -12345.0) score = sfxn.score( xa, xb );
	return d;
}

double
slide_into_contact_and_score_DEPRICATED(
	protocols::sic_dock::SICFast    const & sic,
	protocols::sic_dock::RigidScore const & sfxn,
	core::kinematics::Stub                & xa,
	core::kinematics::Stub          const & xb,
	numeric::xyzVector<core::Real>  const & ori,
	core::Real                            & score
){
	double d = sic.slide_into_contact_DEPRICATED(xa,xb,ori);
	xa.v += d*ori;
	if(score != -12345.0) score = sfxn.score( Xform(xa.M,xa.v), Xform(xb.M,xb.v) );
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

core::Size
count_CBs(
	core::pose::Pose const & pose
){
	platform::Size cbcount = 0;
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


vector1<numeric::xyzVector<core::Real> >
get_CB_Vecs_from_pose(
	core::pose::Pose const & pose
){
	vector1<numeric::xyzVector<core::Real> > CBs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		if( pose.residue(ir).has("CB") ){
			CBs.push_back( pose.residue(ir).xyz("CB") );
		}
	}
	return CBs;
}
vector1<numeric::xyzVector<core::Real> >
get_CB_Vecs_from_map(
	core::pose::Pose const & pose,
	core::id::AtomID_Map<core::Real> const & map
){
	vector1<numeric::xyzVector<core::Real> > CBs;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia) {
			Real wt = fabs(map[AtomID(ia,ir)]);
			if( wt > 0.0001 ) CBs.push_back( pose.residue(ir).xyz(ia) );
		}
	}
	return CBs;
}

vector1<core::Real>
cb_weights_from_pose(
	core::pose::Pose const & pose
){
	vector1<core::Real> wts;
	for(Size i = 1; i <= pose.n_residue(); ++i){
		if(pose.residue(i).has("CB")) {
			wts.push_back( cb_weight(pose,i) );
		}
	}
	return wts;
}
vector1<core::Real>
cb_weights_from_map(
	core::pose::Pose const & pose,
	core::id::AtomID_Map<core::Real> const & map
){
	vector1<core::Real> wts;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		for(Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia) {
			Real wt = fabs(map[AtomID(ia,ir)]);
			if( wt > 0.0001 ) wts.push_back(wt);
		}
	}
	return wts;
}

vector1<core::Size> range(core::Size beg, core::Size end){
	vector1<core::Size> v;
	for(core::Size i = beg; i < end; ++i) v.push_back(i);
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
// 	vector1<Real> rsd_sasa;
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
// 	// core::Real nnt=0.0,nct=0.0,gnt=0.0,gct=0.0;
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

bool residue_is_floppy(core::pose::Pose const & pose, Size const ir, Real const ttrim_cut, Size const nfold){
	if(pose.n_residue()==1) return false;
	Real contacts = 0.0;
	for(Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia){
		Vec const & p1 = pose.residue(ir).xyz(ia);
		for(Size jr = 1; jr <= pose.n_residue(); ++jr){
			if( abs((int)ir-(int)jr) < 3 ) continue;
			for(Size ja = 1; ja <= pose.residue(jr).nheavyatoms(); ++ja){
				Vec const & p2 = pose.residue(jr).xyz(ja);
				if(p1.distance_squared(p2) < 25.0) contacts += 1.0;
			}
		}
		// !!!!!!!!!! assume "standard" bb alignment
		for(Size ifold = 1; ifold < nfold; ++ifold){
			Mat R = rotation_matrix_degrees(Vec(0,0,1),360.0/(Real)ifold/(Real)nfold);
			for(Size jr = 1; jr <= pose.n_residue(); ++jr){
				for(Size ja = 1; ja <= pose.residue(jr).nheavyatoms(); ++ja){
					Vec const & p2 = R * pose.residue(jr).xyz(ja);
					if(p1.distance_squared(p2) < 25.0) contacts += 1.0;
				}
			}
		}
	}
	contacts /= (Real)pose.residue(ir).nheavyatoms();
	return contacts < ttrim_cut;
}

void
auto_trim_floppy_termini(core::pose::Pose & pose, Real const ttrim_cut, Size const nfold){
	Size ntrim=0,ctrim=0;
	while( residue_is_floppy(pose,       1        ,ttrim_cut,nfold) ){
		pose.conformation().delete_residue_slow(       1        );
		++ntrim;
	}
	TR << "trimmed " << ntrim << " n-terminal residues" << endl;
	while( residue_is_floppy(pose,pose.n_residue(),ttrim_cut,nfold) ){
		pose.conformation().delete_residue_slow(pose.n_residue());
		++ctrim;
	}
	TR << "trimmed " << ctrim << " c-terminal residues" << endl;
}

void dump_points_pdb(vector1<Vec> const & p, std::string fn) {
	using namespace ObjexxFCL::format;
	utility::io::ozstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		std::string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}
void dump_points_pdb(vector1<Vec> const & p, Vec t, std::string fn) {
	using namespace ObjexxFCL::format;
	utility::io::ozstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		std::string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x()+t.x())<<F(8,3,p[i].y()+t.y())<<F(8,3,p[i].z()+t.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}
void trans_pose( Pose & pose, Vec const & trans, Size start, Size end) {
	if(0==end) end = pose.n_residue();
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + (Vec)trans );
		}
	}
}
void rot_pose( Pose & pose, Mat const & rot, Size start, Size end ) {
	if(0==end) end = pose.n_residue();
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}
void rot_pose( Pose & pose, Mat const & rot, Vec const & cen, Size start, Size end ) {
	trans_pose(pose,-cen,start,end);
	rot_pose(pose,rot,start,end);
	trans_pose(pose,cen,start,end);
}
void rot_pose( Pose & pose, Vec const & axis, Real const & ang, Size start, Size end ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),start,end);
}
void rot_pose( Pose & pose, Vec const & axis, Real const & ang, Vec const & cen, Size start, Size end) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen,start,end);
}
void alignaxis(Pose & pose, Vec newaxis, Vec oldaxis, Vec cen ) {
	newaxis.normalize();
	oldaxis.normalize();
	if(fabs(newaxis.dot(oldaxis)) < 0.9999) {
		Vec axis = newaxis.cross(oldaxis).normalized();
		Real ang = (Real)-acos(numeric::max((Real)-1.0,numeric::min((Real)1.0,newaxis.dot(oldaxis))))*(Real)180.0/numeric::constants::f::pi;
		rot_pose(pose,axis,ang,cen);
	}
}
numeric::xyzTransform<core::Real> alignaxis_xform (numeric::xyzVector<core::Real>  newaxis, numeric::xyzVector<core::Real>  oldaxis, numeric::xyzVector<core::Real>  cen){
	newaxis.normalize();
	oldaxis.normalize();
	if(fabs(newaxis.dot(oldaxis)) < 0.9999) {
		Vec axis = newaxis.cross(oldaxis).normalized();
		Real ang = (Real)-acos(numeric::max((Real)-1.0,numeric::min((Real)1.0,newaxis.dot(oldaxis))))*(Real)180.0/numeric::constants::f::pi;
		return numeric::xyzTransform<core::Real>::rot_deg(axis,ang,cen);
	}
	return numeric::xyzTransform<core::Real>();
}

void xform_pose( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres, Size eres ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
    }
  }
}
void xform_pose_rev( core::pose::Pose & pose, core::kinematics::Stub const & s ) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
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

Vec center_of_geom(core::pose::Pose const & pose, Size str, Size end) {
	if(!end) end = pose.n_residue();
	Vec c(0,0,0);
	for(Size i = str; i <= end; ++i) {
		c += pose.xyz(AtomID(2,i));
	}
	c /= Real(end-str+1);
	return c;
}




core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::Residue const & fixd_resi) {
  core::kinematics::Stub s;
  s.M = alignVectorSets(move_resi.xyz(1)-move_resi.xyz(2),move_resi.xyz(3)-move_resi.xyz(2),fixd_resi.xyz(1)-fixd_resi.xyz(2),fixd_resi.xyz(3)-fixd_resi.xyz(2));
  s.v = fixd_resi.xyz(2)-s.M*move_resi.xyz(2);
  return s;
}

Real brute_mindis(vector1<Vec> const & pa, vector1<Vec> const & pb, Vec const & ofst=Vec(0,0,0)) {
	Real mindis = 9e9;
	for(vector1<Vec>::const_iterator i = pa.begin(); i != pa.end(); ++i) {
		for(vector1<Vec>::const_iterator j = pb.begin(); j != pb.end(); ++j) {
			mindis = min( mindis, i->distance_squared(*j+ofst) );
		}
	}
	return mindis;
}

} // sic_dock
} // protocols
