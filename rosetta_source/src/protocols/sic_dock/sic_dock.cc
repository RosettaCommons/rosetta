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

	inline double sigmoid( double const & sqdist, double const & start, double const & stop ) {
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

  SICFast::~SICFast() {
    if(xh2_bb_) delete xh2_bb_;
    if(xh2_cb_) delete xh2_cb_;
  }

  SICFast::SICFast() :
		CTD(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]()),
		CLD(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]()),
		CTD2(sqr(CTD)),CLD2(sqr(CLD)),BIN(CLD*basic::options::option[basic::options::OptionKeys::sicdock::hash_2D_vs_3D]()),
		xh2_bb_(NULL),xh2_cb_(NULL)
	{}



	// void SICFast::init(core::pose::Pose const & cmp1in, vector1<Vec> cmp1cbs, vector1<double> const & cmp1wts,
	// 									 core::pose::Pose const & cmp2in, vector1<Vec> cmp2cbs, vector1<double> const & cmp2wts
	// 									 )
	// {
	// 	xh2_bb_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]()+0.5);
	// 	xh2_cb_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
	// 	xh2_bb_->init_with_pose(cmp2in,BB);
	// 	xh2_cb_->init_with_pose(cmp2in,cmp2wts,CB);
	// 	w1_ = cmp1wts;
	// 	w2_ = cmp2wts;
	// }

	void SICFast::init(core::pose::Pose const & pose1,
										 core::pose::Pose const & pose2,
										 core::id::AtomID_Map<core::Real> const & clash_atoms1,
										 core::id::AtomID_Map<core::Real> const & clash_atoms2,
										 core::id::AtomID_Map<core::Real> const & score_atoms1,
										 core::id::AtomID_Map<core::Real> const & score_atoms2
										 )
	{
		using core::id::AtomID;
		xh2_bb_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::clash_dis]()+0.5);
		xh2_cb_ = new xyzStripeHashPose(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
		xh2_bb_->init_with_pose(pose2,clash_atoms2);
		xh2_cb_->init_with_pose(pose2,score_atoms2);
		for(int ir = 1; ir <= pose1.n_residue(); ++ir) {
			for(int ia = 1; ia <= clash_atoms1.n_atom(ir); ++ia) { if(clash_atoms1[AtomID(ia,ir)] > 0) clash1_.push_back(pose1.xyz(AtomID(ia,ir))); }
			for(int ia = 1; ia <= score_atoms1.n_atom(ir); ++ia) { if(score_atoms1[AtomID(ia,ir)] > 0) score1_.push_back(pose1.xyz(AtomID(ia,ir))); }
			for(int ia = 1; ia <= score_atoms1.n_atom(ir); ++ia) { if(score_atoms1[AtomID(ia,ir)] > 0) w1_.push_back(score_atoms1[AtomID(ia,ir)]); }
		}
		for(int ir = 1; ir <= pose2.n_residue(); ++ir) {
			for(int ia = 1; ia <= clash_atoms2.n_atom(ir); ++ia) { if(clash_atoms2[AtomID(ia,ir)] > 0) clash2_.push_back(pose2.xyz(AtomID(ia,ir))); }
			for(int ia = 1; ia <= score_atoms2.n_atom(ir); ++ia) { if(score_atoms2[AtomID(ia,ir)] > 0) score2_.push_back(pose2.xyz(AtomID(ia,ir))); }
			for(int ia = 1; ia <= score_atoms2.n_atom(ir); ++ia) { if(score_atoms2[AtomID(ia,ir)] > 0) w2_.push_back(score_atoms2[AtomID(ia,ir)]); }
		}
	}

	double SICFast::slide_into_contact(core::kinematics::Stub const & xa,
																		 core::kinematics::Stub const & xb,
																		 utility::vector1<Vec>          pa,
																		 utility::vector1<Vec>          pb,
																		 utility::vector1<Vec>  const & cba,
																		 utility::vector1<Vec>  const & cbb,
																		 Vec                            ori,
																		 double                       & score
																		 )
	{
		rotate_points(pa,pb,ori);
		get_bounds(pa,pb);
		fill_plane_hash(pa,pb);
		double const mindis_approx = get_mindis_with_plane_hashes();
		double const mindis = refine_mindis_with_xyzHash(xa,xb,pa,pb,mindis_approx,ori);
		//cerr << brute_mindis(pa,pb,Vec(0,0,-mindis)) << endl;
		// if( fabs(CLD2-brute_mindis(pa,pb,Vec(0,0,-mindis))) > 0.0001 ) utility_exit_with_message("DIAF!");
		if(score != -12345.0) score = get_score(xa,xb,cba,cbb,ori,mindis);
		return mindis;
	}



 	void SICFast::rotate_points(vector1<Vec> & pa, vector1<Vec> & pb, Vec ori) {
		// // get points, rotated ro ori is 0,0,1, might already be done
		// Mat rot = Mat::identity();
		// if     ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), (double)-acos(Vec(0,0,1).dot(ori)) );
		// else if( ori.dot(Vec(0,0,1)) <  0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), (double)-acos(Vec(0,0,1).dot(ori)) );
		// if( rot != Mat::identity() ) {
		//         for(vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
		//         for(vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
		// }
		Mat rot = rotation_matrix_degrees( (ori.z() < -0.99999) ? Vec(1,0,0) : (Vec(0,0,1)+ori.normalized())/2.0 , 180.0 );
		for(vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
		for(vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	}
	void SICFast::get_bounds(vector1<Vec> & pa, vector1<Vec> & pb) {
		// get bounds for plane hashes
		xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
		for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
			xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
			ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
		}
		for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
			xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
			ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
		}
		xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
		ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);
	}
	double SICFast::get_score(core::kinematics::Stub const & xa, core::kinematics::Stub const & xb,
														vector1<Vec> const & cba, vector1<Vec> const & cbb, Vec ori, double mindis)
	{
		vector1<double> const & wa = w2_;
		vector1<double> const & wb = w1_;
		numeric::geometry::hashing::xyzStripeHash<double> const * xh( xh2_cb_);
		//		Mat R = numeric::rotation_matrix_degrees(axsa,-anga);
		float score = 0.0;
		vector1<double>::const_iterator iwb = wb.begin();
		for(vector1<Vec>::const_iterator i = cbb.begin(); i != cbb.end(); ++i,++iwb) {
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
							score += sigmoid(d2, CLD, CTD ) * a2.w * (*iwb);
						}
					}
				}
			}
		}
		// double cbcount = 0.0;
		// vector1<double>::const_iterator iwa = wa.begin();
		// for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia,++iwa) {
		// 	vector1<double>::const_iterator iwb = wb.begin();
		// 	for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib,++iwb) {
		// 		double d2 = ib->distance_squared( (*ia) + (mindis*ori) );
		// 		if( d2 < CTD2 ) {
		// 			cbcount += sigmoid(d2, CLD, CTD ) * (*iwa) * (*iwb);
		// 		}
		// 	}
		// }
		// if( fabs(cbcount-score) > 0.0001 ) {
		// 	cout << "btute/hash score mismatch!!!! " << anga << " " << angb << " " << cbcount << " " << score << endl;
		// }
		return score;
	}
	double SICFast::refine_mindis_with_xyzHash(core::kinematics::Stub const & xa, core::kinematics::Stub const & xb,
																						 vector1<Vec> const & pa, vector1<Vec> const & pb, double mindis, Vec ori)
	{

		numeric::geometry::hashing::xyzStripeHash<double> const * xh( xh2_bb_ );
		Mat Rori = rotation_matrix_degrees( (ori.z() < -0.99999) ? Vec(1,0,0) : (Vec(0,0,1)+ori.normalized())/2.0 , 180.0 );
		//Mat R    = numeric::rotation_matrix_degrees(axsa,-anga);
		//Mat Rinv = numeric::rotation_matrix_degrees(axsa, anga);
		while(true){
			double correction_hash = 9e9;
			for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
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
		// double correction_safe = 9e9;
		// for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
		// 	Vec const v = ((*ib)-Vec(0,0,mindis));
		// 	// dbg2.push_back(v);
		// 	for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
		// 		// correction_safe = min(correction_safe,ia->distance_squared(v));
		// 		double const dxy2 = (ia->x()-v.x())*(ia->x()-v.x()) + (ia->y()-v.y())*(ia->y()-v.y());
		// 		if( dxy2 >= CLD2 ) continue;
		// 		// cout << "BRUTE " << dxy2 << endl;
		// 		double const dz = v.z() - ia->z() - sqrt(CLD2-dxy2);
		// 		if( dz < correction_safe) correction_safe = dz;
		// 	}
		// }
		// #ifdef USE_OPENMP
		// #pragma omp critical
		// #endif
		// if( fabs(correction_safe-correction_hash) > 0.01 ) {
		// 	cout << F(9,5,correction_hash) << " " << F(9,5,correction_safe) << " " << mindis << endl;
		// 	// utility_exit_with_message("FOO");
		// }

		return mindis;
	}

	void SICFast::fill_plane_hash(vector1<Vec> & pa, vector1<Vec> & pb) {
		xlb = (int)(xmn/BIN)-2; xub = (int)(xmx/BIN+0.999999999)+2; // one extra on each side for correctness,
		ylb = (int)(ymn/BIN)-2; yub = (int)(ymx/BIN+0.999999999)+2; // and one extra for outside atoms
		ha.dimension(xub-xlb+1,yub-ylb+1,Vec(0,0,-9e9));
		hb.dimension(xub-xlb+1,yub-ylb+1,Vec(0,0, 9e9));
		int const xsize = xub-xlb+1;
		int const ysize = yub-ylb+1;
		for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
			int const ix = (int)((ia->x()/BIN)-xlb+0.999999999);
			int const iy = (int)((ia->y()/BIN)-ylb+0.999999999);
			if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
			if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
			// bool const test = !( ix < 1 || ix > xsize || iy < 1 || iy > ysize) && ha(ix,iy).z() < ia->z();
			// ha(ix,iy) = test ? *ia : ha(ix,iy);
		}
		for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
			int const ix = (int)((ib->x()/BIN)-xlb+0.999999999);
			int const iy = (int)((ib->y()/BIN)-ylb+0.999999999);
			if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
			if( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
			// bool const test = !( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) && hb(ix,iy).z() > ib->z();
			// hb(ix,iy) = test ? *ib : hb(ix,iy);
		}

	}
	double SICFast::get_mindis_with_plane_hashes() {
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
    if(i>1           ) nmark += flood_fill3D(i-1,j  ,k  ,grid,t);
    if(i<grid.size1()) nmark += flood_fill3D(i+1,j  ,k  ,grid,t);
    if(j>1           ) nmark += flood_fill3D(i  ,j-1,k  ,grid,t);
    if(j<grid.size2()) nmark += flood_fill3D(i  ,j+1,k  ,grid,t);
    if(k>1           ) nmark += flood_fill3D(i  ,j  ,k-1,grid,t);
    if(k<grid.size3()) nmark += flood_fill3D(i  ,j  ,k+1,grid,t);
    return nmark;
  }



	void termini_exposed(core::pose::Pose const & pose, bool & ntgood, bool & ctgood ){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		core::id::AtomID_Map<Real> atom_sasa;
		core::id::AtomID_Map<bool> atom_subset;
		utility::vector1<Real> rsd_sasa;
		core::pose::initialize_atomid_map(atom_subset, pose, false);
		for(int i = 2; i <= pose.n_residue()-1; ++i) {
			for(int ia = 1; ia <= pose.residue(i).nheavyatoms(); ++ia) {
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
