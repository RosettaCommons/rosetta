// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <protocols/sic_dock/sic_dock.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>

#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

namespace protocols {
namespace sic_dock {

using core::Size;
using core::Real;
using std::string;
using utility::vector1;
using ObjexxFCL::fmt::I;
using ObjexxFCL::fmt::F;
using ObjexxFCL::fmt::RJ;
using numeric::min;
using numeric::max;
using std::cout;
using std::cerr;
using std::endl;
typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

template<typename T> inline T sqr(T x) { return x*x; }


void dump_points_pdb(utility::vector1<Vec> const & p, std::string fn) {
        using namespace ObjexxFCL::fmt;
        std::ofstream o(fn.c_str());
        for(Size i = 1; i <= p.size(); ++i) {
                std::string rn = "VIZ";                o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        }        o.close();
}


inline double dist_score( double const & sqdist, double const & start, double const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		double dist = sqrt( sqdist );
		return (stop-dist)/(stop-start);
		//return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
	}
}

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


SICFast::~SICFast() {
	if(xh2c_) delete xh2c_;
	if(xh2s_) delete xh2s_;
}

SICFast::SICFast() :
	CTD(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]()),
	CLD(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]()),
	CTD2(sqr(CTD)),CLD2(sqr(CLD)),BIN(CLD*basic::options::option[basic::options::OptionKeys::sicdock::hash_2D_vs_3D]()),
	xh1c_(NULL),xh1s_(NULL),xh2c_(NULL),xh2s_(NULL)
{}


void 
SICFast::init(
	core::pose::Pose const & pose
) {
	init(pose,pose);
}

void
SICFast::init(
	core::pose::Pose const & pose,
	core::id::AtomID_Map<core::Real> const & clash_atoms,
	core::id::AtomID_Map<core::Real> const & score_atoms
){
	init(pose,pose,clash_atoms,clash_atoms,score_atoms,score_atoms);
}

void 
SICFast::init(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2
){
	using core::id::AtomID;
	core::id::AtomID_Map<core::Real> clashmap1,contactmap1,clashmap2,contactmap2;
	core::pose::initialize_atomid_map(  clashmap1,pose1,-1.0);
	core::pose::initialize_atomid_map(  clashmap2,pose2,-1.0);
	core::pose::initialize_atomid_map(contactmap1,pose1,-1.0);
	core::pose::initialize_atomid_map(contactmap2,pose2,-1.0);
	for(Size i = 1; i <= pose1.n_residue(); ++i) {
		if(pose1.residue(i).has("CB")) 
			contactmap1[AtomID(pose1.residue(i).atom_index("CB"),i)] = min(1.0,(double)neighbor_count(pose1,i)/20.0);
		for(int j = 1; j <= ((pose1.residue(i).name3()=="GLY")?4:5); ++j)
			clashmap1[AtomID(j,i)] = pose1.residue(i).atom_type(j).lj_radius();
	}
	for(Size i = 1; i <= pose2.n_residue(); ++i) {
		if(pose2.residue(i).has("CB"))
			contactmap2[AtomID(pose2.residue(i).atom_index("CB"),i)] = min(1.0,(double)neighbor_count(pose2,i)/20.0);
		for(int j = 1; j <= ((pose2.residue(i).name3()=="GLY")?4:5); ++j)
			clashmap2[AtomID(j,i)] = pose2.residue(i).atom_type(j).lj_radius();
	}
	init(pose1,pose2,clashmap1,clashmap2,contactmap1,contactmap2);
}

void
SICFast::init(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::id::AtomID_Map<core::Real> const & clash_atoms1,
	core::id::AtomID_Map<core::Real> const & clash_atoms2,
	core::id::AtomID_Map<core::Real> const & score_atoms1,
	core::id::AtomID_Map<core::Real> const & score_atoms2
){
	using core::id::AtomID;
	xh1c_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]()+0.5);
	xh2c_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]()+0.5);
	xh1s_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
	xh2s_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
	xh1c_->init_with_pose(pose1,clash_atoms1);
	xh2c_->init_with_pose(pose2,clash_atoms2);
	xh1s_->init_with_pose(pose1,score_atoms1);
	xh2s_->init_with_pose(pose2,score_atoms2);
	for(int i=0;i<xh1s_->natom();++i) w1_.push_back( xh1s_->grid_atoms()[i].w );
	for(int i=0;i<xh2s_->natom();++i) w2_.push_back( xh2s_->grid_atoms()[i].w );
}

double
SICFast::slide_into_contact(
	core::kinematics::Stub const & xa,
	core::kinematics::Stub const & xb,
	Vec                            ori,
	double                       & score
){
	utility::vector1<Vec> pb,pa,sb,sa;
	for(int i=0;i<xh1c_->natom();++i)pa.push_back(xb.local2global(Vec(xh1c_->grid_atoms()[i].x,xh1c_->grid_atoms()[i].y,xh1c_->grid_atoms()[i].z)-xh1c_->translation()));
	for(int i=0;i<xh2c_->natom();++i)pb.push_back(xa.local2global(Vec(xh2c_->grid_atoms()[i].x,xh2c_->grid_atoms()[i].y,xh2c_->grid_atoms()[i].z)-xh2c_->translation()));
	for(int i=0;i<xh1s_->natom();++i)sa.push_back(xb.local2global(Vec(xh1s_->grid_atoms()[i].x,xh1s_->grid_atoms()[i].y,xh1s_->grid_atoms()[i].z)-xh1s_->translation()));
	for(int i=0;i<xh2s_->natom();++i)sb.push_back(xa.local2global(Vec(xh2s_->grid_atoms()[i].x,xh2s_->grid_atoms()[i].y,xh2s_->grid_atoms()[i].z)-xh2s_->translation()));
	rotate_points(pb,pa,ori);
	get_bounds(pb,pa);
	fill_plane_hash(pb,pa);
	double const mindis_approx = get_mindis_with_plane_hashes();
	double const mindis = refine_mindis_with_xyzHash(xa,xb,pb,pa,mindis_approx,ori);
	//cerr << brute_mindis(pb,pa,Vec(0,0,-mindis)) << endl;
	// if( fabs(CLD2-brute_mindis(pb,pa,Vec(0,0,-mindis))) > 0.0001 ) utility_exit_with_message("DIAF!");
	if(score != -12345.0) score = get_score(xa,xb,sb,sa,ori,mindis);
	return mindis;
}



void
SICFast::rotate_points(
	vector1<Vec> & pb,
	vector1<Vec> & pa,
	Vec ori
){
	Mat rot = rotation_matrix_degrees( (ori.z() < -0.99999) ? Vec(1,0,0) : (Vec(0,0,1)+ori.normalized())/2.0 , 180.0 );
	for(vector1<Vec>::iterator ia = pb.begin(); ia != pb.end(); ++ia) *ia = rot*(*ia);
	for(vector1<Vec>::iterator ib = pa.begin(); ib != pa.end(); ++ib) *ib = rot*(*ib);
}

void
SICFast::get_bounds(
	vector1<Vec> const & pb,
	vector1<Vec> const & pa
){
	// get bounds for plane hashes
	xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
	for(vector1<Vec>::const_iterator ia = pb.begin(); ia != pb.end(); ++ia) {
		xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
		ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
	}
	for(vector1<Vec>::const_iterator ib = pa.begin(); ib != pa.end(); ++ib) {
		xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
		ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
	}
	xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
	ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);
}


double
SICFast::get_score(
	core::kinematics::Stub const & xa,
	core::kinematics::Stub const & /*xb*/,
	vector1<Vec> const & /*sb*/,
	vector1<Vec> const & sa,
	Vec ori,
	double mindis
){
	vector1<double> const & wa = w2_;
	vector1<double> const & wb = w1_;
	numeric::geometry::hashing::xyzStripeHash<double> const * xh( xh2s_);
	//		Mat R = numeric::rotation_matrix_degrees(axsa,-anga);
	float score = 0.0;
	vector1<double>::const_iterator iwb = wb.begin();
	for(vector1<Vec>::const_iterator i = sa.begin(); i != sa.end(); ++i,++iwb) {
		Vec v = xa.global2local((*i)-mindis*ori) + xh->translation();
		if( v.x() < -xh->grid_size_ || v.y() < -xh->grid_size_ || v.z() < -xh->grid_size_ ) continue; // worth it?
		if( v.x() >  xh->xmx_       || v.y() >  xh->ymx_       || v.z() >  xh->zmx_       ) continue; // worth it?
		int const ix	= (v.x()<0.0) ? 0 : (int)(numeric::min(xh->xdim_-1,(int)(v.x()/xh->grid_size_)));
		int const iy0	= (v.y()<0.0) ? 0 : (int)(v.y()/xh->grid_size_);
		int const iz0	= (v.z()<0.0) ? 0 : (int)(v.z()/xh->grid_size_);
		int const iyl = numeric::max(0,iy0-1);
		int const izl = numeric::max(0,iz0-1);
		int const iyu = numeric::min((int)xh->ydim_,     iy0+2);
		int const izu = numeric::min((int)xh->zdim_,(int)iz0+2);
		for(int iy = iyl; iy < iyu; ++iy) {
			for(int iz = izl; iz < izu; ++iz) {
				int const ig = ix+xh->xdim_*iy+xh->xdim_*xh->ydim_*iz;
				assert(ig < xh->xdim_*xh->ydim_*xh->zdim_ && ix < xh->xdim_ && iy < xh->ydim_ && iz < xh->zdim_);
				int const igl = xh->grid_stripe_[ig].x;
				int const igu = xh->grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					numeric::geometry::hashing::xyzStripeHash<double>::float4 const & a2 = xh->grid_atoms_[i];
					float const d2 = (v.x()-a2.x)*(v.x()-a2.x) + (v.y()-a2.y)*(v.y()-a2.y) + (v.z()-a2.z)*(v.z()-a2.z);
					if( d2 <= xh->grid_size2_ ) {
						score += dist_score(d2, CLD, CTD ) * a2.w * (*iwb);
					}
				}
			}
		}
	}
	return score;
}

double
SICFast::refine_mindis_with_xyzHash(
	core::kinematics::Stub const & xa,
	core::kinematics::Stub const & /*xb*/,
	vector1<Vec> const & /*pb*/,
	vector1<Vec> const & pa,
	double mindis,
	Vec ori
){
	numeric::geometry::hashing::xyzStripeHash<double> const * xh( xh2c_ );
	Mat Rori = rotation_matrix_degrees( (ori.z() < -0.99999) ? Vec(1,0,0) : (Vec(0,0,1)+ori.normalized())/2.0 , 180.0 );
	//Mat R    = numeric::rotation_matrix_degrees(axsa,-anga);
	//Mat Rinv = numeric::rotation_matrix_degrees(axsa, anga);
	while(true){
		double correction_hash = 9e9;
		for(vector1<Vec>::const_iterator ib = pa.begin(); ib != pa.end(); ++ib) {
			Vec const v = xa.global2local(Rori*((*ib)-Vec(0,0,mindis))) + xh->translation();
			Vec const b = Rori * xa.local2global(v);
			if( v.x() < -xh->grid_size_ || v.y() < -xh->grid_size_ || v.z() < -xh->grid_size_ ) continue; // worth it?
			if( v.x() >  xh->xmx_       || v.y() >  xh->ymx_       || v.z() >  xh->zmx_       ) continue; // worth it?
			int const ix  = (v.x()<0) ? 0 : numeric::min(xh->xdim_-1,(int)(v.x()/xh->grid_size_));
			int const iy0 = (int)((v.y()<0) ? 0 : v.y()/xh->grid_size_);
			int const iz0 = (int)((v.z()<0) ? 0 : v.z()/xh->grid_size_);
			int const iyl = numeric::max(0,iy0-1);
			int const izl = numeric::max(0,iz0-1);
			int const iyu = numeric::min((int)xh->ydim_,     iy0+2);
			int const izu = numeric::min((int)xh->zdim_,(int)iz0+2);
			for(int iy = iyl; iy < iyu; ++iy) {
				for(int iz = izl; iz < izu; ++iz) {
					int const ig = ix+xh->xdim_*iy+xh->xdim_*xh->ydim_*iz;
					assert(ig < xh->xdim_*xh->ydim_*xh->zdim_ && ix < xh->xdim_ && iy < xh->ydim_ && iz < xh->zdim_);
					int const igl = xh->grid_stripe_[ig].x;
					int const igu = xh->grid_stripe_[ig].y;
					for(int i = igl; i < igu; ++i) {
						numeric::geometry::hashing::xyzStripeHash<double>::float4 const & a2 = xh->grid_atoms_[i];
						float const d2 = (v.x()-a2.x)*(v.x()-a2.x) + (v.y()-a2.y)*(v.y()-a2.y) + (v.z()-a2.z)*(v.z()-a2.z);
						Vec const a = Rori * xa.local2global(Vec(a2.x,a2.y,a2.z));
						double const dxy2 = (a.x()-b.x())*(a.x()-b.x()) + (a.y()-b.y())*(a.y()-b.y());
						if( dxy2 >= CLD2 ) continue;
						double const dz = b.z() - a.z() - sqrt(CLD2-dxy2);
						// cout << "HASH " << dz << endl;
						correction_hash = min(dz,correction_hash);
					}
				}
			}
		}
		mindis += correction_hash;
		if( fabs(correction_hash) < 0.001 ) break;
	}
	return mindis;
}

void
SICFast::fill_plane_hash(
	vector1<Vec> const & pb,
	vector1<Vec> const & pa
){
	xlb = (int)(xmn/BIN)-2; xub = (int)(xmx/BIN+0.999999999)+2; // one extra on each side for correctness,
	ylb = (int)(ymn/BIN)-2; yub = (int)(ymx/BIN+0.999999999)+2; // and one extra for outside atoms
	ha.dimension(xub-xlb+1,yub-ylb+1,Vec(0,0,-9e9));
	hb.dimension(xub-xlb+1,yub-ylb+1,Vec(0,0, 9e9));
	int const xsize = xub-xlb+1;
	int const ysize = yub-ylb+1;
	for(vector1<Vec>::const_iterator ia = pb.begin(); ia != pb.end(); ++ia) {
		int const ix = (int)((ia->x()/BIN)-xlb+0.999999999);
		int const iy = (int)((ia->y()/BIN)-ylb+0.999999999);
		if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
		// bool const test = !( ix < 1 || ix > xsize || iy < 1 || iy > ysize) && ha(ix,iy).z() < ia->z();
		// ha(ix,iy) = test ? *ia : ha(ix,iy);
	}
	for(vector1<Vec>::const_iterator ib = pa.begin(); ib != pa.end(); ++ib) {
		int const ix = (int)((ib->x()/BIN)-xlb+0.999999999);
		int const iy = (int)((ib->y()/BIN)-ylb+0.999999999);
		if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
		// bool const test = !( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) && hb(ix,iy).z() > ib->z();
		// hb(ix,iy) = test ? *ib : hb(ix,iy);
	}

}

double
SICFast::get_mindis_with_plane_hashes(
){
	int const xsize=xub-xlb+1, ysize=yub-ylb+1;
	int imna=0,jmna=0,imnb=0,jmnb=0;
	double m = 9e9;
	for(int i = 1; i <= xsize; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
		for(int j = 1; j <= ysize; ++j) {
			for(int k = -1; k <= 1; ++k) {
				if(i+k < 1 || i+k > xsize) continue;
				for(int l = -1; l <= 1; ++l) {
					if(j+l < 1 || j+l > ysize) continue;
					double const xa1=ha(i,j).x(),ya1=ha(i,j).y(),xb1=hb(i+k,j+l).x(),yb1=hb(i+k,j+l).y(),d21=(xa1-xb1)*(xa1-xb1)+(ya1-yb1)*(ya1-yb1);
					if(d21<CLD2){ double const dz=hb(i+k,j+l).z()-ha(i,j).z()-sqrt(CLD2-d21); if(dz<m) m=dz; }
				}
			}
		}
	}
	return m;
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



void
termini_exposed(
	core::pose::Pose const & pose,
	bool & ntgood,
	bool & ctgood
){
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_subset;
	utility::vector1<Real> rsd_sasa;
	core::pose::initialize_atomid_map(atom_subset, pose, false);
	for(int i = 2; i <= (int)pose.n_residue()-1; ++i) {
		for(int ia = 1; ia <= (int)pose.residue(i).nheavyatoms(); ++ia) {
			if(pose.residue(i).atom_is_backbone(ia))
				atom_subset[core::id::AtomID(ia,i)] = true;
		}
	}
	atom_subset[core::id::AtomID(1,1)] = true;
	atom_subset[core::id::AtomID(3,pose.n_residue())] = true;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa,rsd_sasa, 4.0, false, atom_subset );
	Real nexpose = atom_sasa[core::id::AtomID(1,        1       )] / 12.56637 / 5.44 / 5.44;
	Real cexpose = atom_sasa[core::id::AtomID(3,pose.n_residue())] / 12.56637 / 5.44 / 5.44;

	Vec nt = pose.residue(        1       ).xyz("N");
	Vec ct = pose.residue(pose.n_residue()).xyz("C");
	Real nang = angle_degrees(nt,Vec(0,0,0),Vec(nt.x(),nt.y(),0));
	Real cang = angle_degrees(ct,Vec(0,0,0),Vec(ct.x(),ct.y(),0));
	ntgood = nexpose > option[sicdock::term_min_expose]() && nang < option[sicdock::term_max_angle]();
	ctgood = cexpose > option[sicdock::term_min_expose]() && cang < option[sicdock::term_max_angle]();
	// core::Real nnt=0.0,nct=0.0,gnt=0.0,gct=0.0;
	// for(int ir=1; ir<=pose.n_residue(); ++ir) {
	// 	for(int ia=1; ia<=5; ++ia) {
	// 		Vec x = pose.residue(ir).xyz(ia);
	// 		if(angle_degrees(x,Vec(0,0,0),nt) < 15.0 &&  ) {
	// 			nnt += 1.0;
	// 			if( nt.normalized().dot(x) < nt.length() )
	// 				gnt += 1.0;
	// 		}
	// 		if(angle_degrees(x,Vec(0,0,0),ct) < 15.0 ) {
	// 			nct += 1.0;
	// 			if( ct.normalized().dot(x) < ct.length() )
	// 				gct += 1.0;
	// 		}
	// 	}
	// }
}

} // namespace sic_dock
} // namespace protocols
