/// @file
/// @brief


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <math.h>

#include <sstream>
#include <string>
#include <queue>


using namespace basic;
using namespace core;
using namespace core::pose;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

static basic::Tracer TR("cryst.gen");

////////////////////////////////////////////////
// helper functions
inline int pos_mod(int x,int y) {
	int r=x%y; if (r<0) r+=y;
	return r;
}
inline Real pos_mod(Real x,Real y) {
	Real r=std::fmod(x,y); if (r<0) r+=y;
	return r;
}
inline int min_mod(int x,int y) {
	int r=x%y; if (r<-y/2) r+=y;if (r>=y/2) r-=y;
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y); if (r<-0.5*y) r+=y;if (r>=0.5*y) r-=y;
	return r;
}

inline double sign(double x) {return (x >= 0.0) ? 1.0 : -1.0;}
inline double norm4(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

inline bool transforms_equiv_mod1(
		numeric::xyzMatrix<Real> const &S1, numeric::xyzVector<Real> const &T1,
		numeric::xyzMatrix<Real> const &S2, numeric::xyzVector<Real> const &T2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() ) +
		std::fabs( min_mod( T1[0] - T2[0] , 1.0)) +
		std::fabs( min_mod( T1[1] - T2[1] , 1.0))  +
		std::fabs( min_mod( T1[2] - T2[2] , 1.0));
	return (err <= 1e-6);
}
inline bool transforms_equiv(
		numeric::xyzMatrix<Real> const &S1, numeric::xyzVector<Real> const &T1,
		numeric::xyzMatrix<Real> const &S2, numeric::xyzVector<Real> const &T2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() ) +
		std::fabs( ( T1[0] - T2[0] )) +
		std::fabs( ( T1[1] - T2[1] ))  +
		std::fabs( ( T1[2] - T2[2] ));
	return (err <= 1e-6);
}
inline bool transforms_equiv(
		numeric::xyzMatrix<Real> const &S1,
		numeric::xyzMatrix<Real> const &S2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() );
	return (err <= 1e-6);
}

void
disp_symmop( core::kinematics::RT const &rt ) {
	numeric::xyzMatrix<Real> const &S=rt.get_rotation();
	numeric::xyzVector<Real> const &T=rt.get_translation();

	std::cerr << "S=["<<S.xx()<<","<<S.xy()<<","<<S.xz()<<";"
		<<S.yx()<<","<<S.yy()<<","<<S.yz()<<";"
		<<S.zx()<<","<<S.zy()<<","<<S.zz()<<"]  T=["
		<<T[0]<<","<<T[1]<<","<<T[2]<<"]" << std::endl;
}



///////////////////////////////////////////////////////////////////////////////
struct Quat {
	Real x_,y_,z_,w_;
	Quat() {
		x_=0; y_=0; z_=0; w_=1;
	}
	Quat( Real x_in,Real y_in,Real z_in,Real w_in) {
		x_=x_in; y_=y_in; z_=z_in; w_=w_in;
	}

	// construct from a rotation matrix
	// should be stable
	Quat( numeric::xyzMatrix<Real> const& R) {
		Real S;
		if (R.xx() + R.yy() + R.zz() > 0) {
			S = sqrt(1 + R.xx() + R.yy() + R.zz()) * 2;
			w_ = 0.25 * S;
			x_ = ( R.zy() - R.yz() ) * S;
			y_ = ( R.xz() - R.zx() ) * S;
			z_ = ( R.yx() - R.xy() ) * S;
		} else if ( R.xx() > R.yy() && R.xx() > R.zz() )  {
			S  = sqrt( 1.0 + R.xx() - R.yy() - R.zz() ) * 2;
			x_ = 0.25 * S;
			y_ = (R.yx() + R.xy() ) / S;
			z_ = (R.zx() + R.xz() ) / S;
			w_ = (R.zy() - R.yz() ) / S;
		} else if ( R.yy() > R.zz() ) {
			S  = sqrt( 1.0 - R.xx() + R.yy() - R.zz() ) * 2;
			x_ = (R.yx() + R.xy() ) / S;
			y_ = 0.25 * S;
			z_ = (R.zy() + R.yz() ) / S;
			w_ = (R.xz() - R.zx() ) / S;
		} else {
			S  = sqrt( 1.0 - R.xx() - R.yy() + R.zz() ) * 2;
			x_ = (R.xz() + R.zx() ) / S;
			y_ = (R.zy() + R.yz() ) / S;
			z_ = 0.25 * S;
			w_ = (R.yx() - R.xy() ) / S;
		}
	}

	numeric::xyzMatrix<Real> asR() const {
		if (1-w_*w_ < 1e-6)
			return numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

		Real xx = x_*x_, xy = x_*y_, xz = x_*z_, xw = x_*w_;
		Real yy = y_*y_, yz = y_*z_, yw = y_*w_;
		Real zz = z_*z_, zw = z_*w_;
		//Real ww = w_*w_;

		return numeric::xyzMatrix<core::Real>::rows(
			1 - 2 * ( yy+zz ) ,     2 * ( xy-zw ) ,     2 * ( xz+yw ) ,
			    2 * ( xy+zw ) , 1 - 2 * ( xx+zz ) ,     2 * ( yz-xw ) ,
			    2 * ( xz-yw ) ,     2 * ( yz+xw ) , 1 - 2 * ( xx+yy ) );
	}
};

struct pointGroupHit {
	pointGroupHit() {
		symmop1=symmop2=symmop3=0;
	}

	pointGroupHit(
			std::string name_in, numeric::xyzVector<Real> axis_in, numeric::xyzVector<Real> axis2_in,
			numeric::xyzVector<Real> origin_in,
			int symmop1_in,
			int symmop2_in  , numeric::xyzVector<Real> ls2_in,
			int symmop3_in=0, numeric::xyzVector<Real> ls3_in=numeric::xyzVector<Real>(0,0,0)
	) {
		name=name_in;
		symmop1=symmop1_in;
		symmop2=symmop2_in;
		symmop3=symmop3_in;
		axis=axis_in;
		axis2=axis2_in;
		ls2=ls2_in;
		origin=origin_in;
		ls3=ls3_in;
	}
	pointGroupHit(
			std::string name_in, numeric::xyzVector<Real> axis_in, numeric::xyzVector<Real> origin_in,
			int symmop1_in,
			int symmop2_in  , numeric::xyzVector<Real> ls2_in
	) {
		name=name_in;
		symmop1=symmop1_in;
		symmop2=symmop2_in;
		symmop3=0;
		axis=axis_in;
		axis2=numeric::xyzVector<Real>(0,0,0);
		ls2=ls2_in;
		origin=origin_in;
		ls3=numeric::xyzVector<Real>(0,0,0);
	}

	std::string name;
	int symmop1, symmop2,symmop3;
	numeric::xyzVector<Real> axis, axis2, origin, ls2, ls3;
	Real radius;
};

struct latticeHit {
	latticeHit( pointGroupHit hit1_in, pointGroupHit hit2_in, Real angle_in, Real offset_in, Real shift_in ) {
		hit1=hit1_in; hit2=hit2_in; angle=angle_in; offset=offset_in; shift=shift_in;
	}

	pointGroupHit hit1,hit2;
	Real angle;
	Real offset;
	Real shift;
};

// get point symm from symmops
utility::vector1< pointGroupHit >
get_point_groups(utility::vector1<core::kinematics::RT> const &rts, numeric::xyzMatrix<core::Real> skewM, bool fastmode ) {
	utility::vector1< pointGroupHit > hits, hits_dedup;
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	int nxforms = (int)rts.size();

	// [1] point symm around (0,0,0)
	int MAXS = 1;
	for (int x=-MAXS; x<=MAXS; ++x)
	for (int y=-MAXS; y<=MAXS; ++y)
	for (int z=-MAXS; z<=MAXS; ++z) {
		for (int i=2; i<=nxforms; ++i) {
			numeric::xyzMatrix<Real> const &Si=rts[i].get_rotation();
			numeric::xyzVector<Real> lattice_shift = numeric::xyzVector<Real>(x,y,z);
			numeric::xyzVector<Real> Ti = rts[i].get_translation() + lattice_shift;

			if (transforms_equiv(Si,identity)) continue;

			Quat Q(skewM*Si*numeric::inverse(skewM));
			numeric::xyzVector<Real> axis(Q.x_, Q.y_, Q.z_);
			axis.normalize();

			// 2/3/4/6
			//int symm = 1;
			numeric::xyzMatrix<Real> ST = Si*Si;
			numeric::xyzVector<Real> TT = Ti + Si*Ti;
			numeric::xyzVector<Real> TR = Ti + TT;
			if (transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0))) {
				hits.push_back( pointGroupHit("C2", axis, TR/2, 1, i, lattice_shift) );
				continue;
			}

			ST = Si*ST;
			TT = Ti+Si*TT;
			TR = TR+TT;
			if (transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0))) {
				hits.push_back( pointGroupHit("C3", axis, TR/3, 1, i, lattice_shift) );
				continue;
			}

			ST = Si*ST;
			TT = Ti+Si*TT;
			TR = TR+TT;
			if (transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0))) {
				hits.push_back( pointGroupHit("C4", axis, TR/4, 1, i, lattice_shift) );
				continue;
			}

			ST = Si*ST;
			TT = Ti+Si*TT;
			TR = TR+TT;
			// dont care about 5-folds

			ST = Si*ST;
			TT = Ti+Si*TT;
			TR = TR+TT;
			if (transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0))) {
				hits.push_back( pointGroupHit("C6", axis, TR/6, 1, i, lattice_shift) );
				continue;
			}
		}
	}

	// expand C2 and Cn groups that are perpendicular
	int nhits = (int)hits.size();
	for (int i=1; i<=nhits; ++i) {
		for (int j=i+1; j<=nhits; ++j) {
			pointGroupHit const &hit_i = hits[i];
			pointGroupHit const &hit_j = hits[j];

			Real dotprod = hit_i.axis.dot(hit_j.axis );
			if ( std::fabs(dotprod) >1e-4) continue;

			if (hit_i.name != "C2" && hit_j.name != "C2") continue;

			numeric::xyzMatrix<Real> const &Si = rts[hit_i.symmop1].get_rotation();
			numeric::xyzMatrix<Real> const &Sj = rts[hit_j.symmop2].get_rotation();
			numeric::xyzVector<Real> const &Ti = rts[hit_i.symmop1].get_translation() + hit_i.ls2;
			numeric::xyzVector<Real> const &Tj = rts[hit_j.symmop2].get_translation() + hit_j.ls2;

			numeric::xyzMatrix<Real> Sij = Si*Sj;
			numeric::xyzVector<Real> Tij = Si*Tj + Ti;
			numeric::xyzMatrix<Real> Sji = Sj*Si;
			numeric::xyzVector<Real> Tji = Sj*Ti + Tj;

			if (hit_i.origin != hit_j.origin) continue;

			pointGroupHit hit_new;
			if (hit_i.name == "C2" && hit_j.name == "C2") {
				hits.push_back( pointGroupHit("D2", hit_i.axis, hit_j.axis, hit_i.origin , 1,hit_i.symmop2, hit_i.ls2, hit_j.symmop2, hit_j.ls2 ) );
			} else if (hit_i.name == "C2") {
				std::string symmgp = "D"+hit_j.name.substr(1,1);
				hits.push_back( pointGroupHit(symmgp, hit_j.axis, hit_i.axis,  hit_i.origin, 1,hit_i.symmop2, hit_i.ls2, hit_j.symmop2, hit_j.ls2) );
			} else if (hit_j.name == "C2") {
				std::string symmgp = "D"+hit_i.name.substr(1,1);
				hits.push_back( pointGroupHit(symmgp, hit_i.axis, hit_j.axis, hit_i.origin, 1,hit_i.symmop2, hit_i.ls2, hit_j.symmop2, hit_j.ls2) );
			}
		}
	}

	// find T's and O's
	for (int i=1; i<=nhits; ++i) {
		for (int j=i+1; j<=nhits; ++j) {
			pointGroupHit const &hit_i = hits[i];
			pointGroupHit const &hit_j = hits[j];

			if (hit_i.name != "C3" && hit_j.name != "C3") continue;

			Real b = hit_i.axis.dot(hit_j.axis );
			if ( std::fabs(b) > 0.99999) continue;

			numeric::xyzVector<Real> w0=hit_i.origin-hit_j.origin;
			Real d=dot(hit_j.axis,w0);
			Real e=dot(hit_i.axis,w0);
			numeric::xyzVector<Real> x = (hit_i.origin-hit_j.origin)+((b*e-d)*hit_i.axis-(e-b*d)*hit_j.axis)/(1-b*b);
			Real offset = x.length();
			Real sc=(b*e-d)/(1-b*b);
			numeric::xyzVector<Real> new_origin=hit_i.origin+hit_i.axis*sc;

			if (offset<1e-4) continue;

			if (std::fabs( b - 1.0/3.0) <= 1e-4) {
				numeric::xyzMatrix<Real> const &Si = rts[hit_i.symmop1].get_rotation();
				numeric::xyzMatrix<Real> const &Sj = rts[hit_j.symmop2].get_rotation();
				numeric::xyzVector<Real> const &Ti = rts[hit_i.symmop1].get_translation() + hit_i.ls2;
				numeric::xyzVector<Real> const &Tj = rts[hit_j.symmop2].get_translation() + hit_j.ls2;

				numeric::xyzMatrix<Real> Sij = Si*numeric::inverse(Sj);
				numeric::xyzVector<Real> Tij = Ti-Tj;
				numeric::xyzMatrix<Real> ST = Sij*Sij;
				numeric::xyzVector<Real> TT = Tij + Sij*Tij;
				numeric::xyzVector<Real> TR = Ti + TT;
				ST = Sij*ST;
				TT = Tij+Sij*TT;
				TR = TR+TT;
				if (transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0))) {
					hits.push_back( pointGroupHit("T", hit_i.axis, hit_j.axis, new_origin, 1,hit_i.symmop2, hit_i.ls2, hit_j.symmop2, hit_j.ls2 ) );
				}
			} else if (std::fabs( b ) <= 1e-4) {
				numeric::xyzMatrix<Real> const &Si = rts[hit_i.symmop1].get_rotation();
				numeric::xyzMatrix<Real> const &Sj = rts[hit_j.symmop2].get_rotation();
				numeric::xyzVector<Real> const &Ti = rts[hit_i.symmop1].get_translation() + hit_i.ls2;
				numeric::xyzVector<Real> const &Tj = rts[hit_j.symmop2].get_translation() + hit_j.ls2;

				numeric::xyzMatrix<Real> Sij = Si*numeric::inverse(Sj);
				numeric::xyzVector<Real> Tij = Ti-Tj;
				numeric::xyzMatrix<Real> ST = Sij*Sij;
				numeric::xyzVector<Real> TT = Tij + Sij*Tij;
				numeric::xyzVector<Real> TR = Ti + TT;
				ST = Sij*ST;
				TT = Tij+Sij*TT;
				TR = TR+TT;
				ST = Sij*ST;
				TT = Tij+Sij*TT;
				TR = TR+TT;
				if (transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0))) {
					hits.push_back( pointGroupHit("O", hit_i.axis, hit_j.axis, new_origin, 1,hit_i.symmop2, hit_i.ls2, hit_j.symmop2, hit_j.ls2 ) );
				}
			}
		}
	}


	// remove duplicates
	nhits = (int)hits.size();
	for (int i=1; i<=nhits; ++i) {
		bool dup=false;
		for (int j=1; j<=(int)hits_dedup.size() && !dup; ++j) {
			Real delaxis = (hits[i].axis - hits_dedup[j].axis).length();
			Real delorigin = (hits[i].origin - hits_dedup[j].origin).length();

			if (hits[i].name == hits_dedup[j].name
					&& delaxis <= 1e-4
					&& delorigin <= 1e-4)
				dup=true;

			if (hits[i].name == "O" && hits_dedup[j].name == "O"
					&& delorigin <=1e-4
					&& std::fabs( hits[i].axis.dot(hits_dedup[j].axis ) - 1.0/3.0) <1e-4 )
				dup=true;

			if (hits[i].name == "T" && hits_dedup[j].name == "T"
					&& delorigin <=1e-4
					&& std::fabs( hits[i].axis.dot(hits_dedup[j].axis )) <1e-4 )
				dup=true;
		}
		if (!dup) hits_dedup.push_back(hits[i]);
	}

	return hits_dedup;
}


bool
check_if_forms_lattice(
		utility::vector1<core::kinematics::RT> const &rts,
		utility::vector1<core::kinematics::RT> const &rts_all,
		bool verbose=false ) {
	int EXPAND_ROUNDS=12;
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);

	// don't include duplicate symmops!
 	utility::vector1<core::kinematics::RT> rts_exp = rts;
 	int nxformsOrig = (int)rts_exp.size();

	for (int rd=1; rd<=EXPAND_ROUNDS; ++rd) {
		int nxforms = (int)rts_exp.size();
		for (int i=1; i<=nxforms; ++i) {
			for (int j=1; j<=nxformsOrig; ++j) {
				numeric::xyzMatrix<Real> const &Si=rts_exp[i].get_rotation();
				numeric::xyzVector<Real> const &Ti=rts_exp[i].get_translation();
				numeric::xyzMatrix<Real> const &Sj=rts_exp[j].get_rotation();
				numeric::xyzVector<Real> const &Tj=rts_exp[j].get_translation();

				numeric::xyzMatrix<Real> Sij = Sj*Si;
				numeric::xyzVector<Real> Tij = Sj*Ti + Tj;

				if (std::fabs(Tij[0])>2 || std::fabs(Tij[1])>2 || std::fabs(Tij[2])>2) continue;

				bool unique_xform = true;
				for (int k=1; k<=(int)rts_exp.size(); ++k) {
					if ( transforms_equiv( Sij, Tij, rts_exp[k].get_rotation(), rts_exp[k].get_translation() )) {
						unique_xform = false;
						break;
					}
				}
				if (unique_xform)
					rts_exp.push_back( core::kinematics::RT( Sij, Tij ) );
			}
		}
	}

	//std::cerr << "-----" << std::endl;
	//for (int i=1; i<=(int)rts_exp.size(); ++i) {
	//	numeric::xyzMatrix<Real> const &Si=rts_exp[i].get_rotation();
	//	if (Si == identity)
	//		disp_symmop(rts_exp[i]);
	//}

	// check if:
	//    1 - all generated transformations are in the lattice group
	//    2 - all lattice transforms are generated
	bool connected = true;
	for (int i=1; i<=(int)rts_all.size() && connected; ++i) {
		bool contains_j = false;
		for (int j=1; j<=(int)rts_exp.size(); ++j) {
			bool i_equals_j = transforms_equiv_mod1(
					rts_exp[j].get_rotation(), rts_exp[j].get_translation(),
					rts_all[i].get_rotation(), rts_all[i].get_translation());
			contains_j |= i_equals_j;
		}
		connected &= contains_j;
		//if (!contains_j) std::cerr << "fail to connect " << i << std::endl;
	}
	if(!connected) return false;

	bool has_00p=false, has_00m=false, has_0m0=false, has_0p0=false, has_p00=false, has_m00=false;
	for (int i=1; i<=(int)rts_exp.size(); ++i) {
		numeric::xyzMatrix<Real> const &Si=rts_exp[i].get_rotation();
		numeric::xyzVector<Real> const &Ti=rts_exp[i].get_translation();
		if (transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0, 1))) has_00p = true;
		if (transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0,-1))) has_00m = true;
		if (transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0, 1,0))) has_0p0 = true;
		if (transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,-1,0))) has_0m0 = true;
		if (transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>( 1,0,0))) has_p00 = true;
		if (transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(-1,0,0))) has_m00 = true;
	}

	connected &= has_00p;
	connected &= has_00m;
	connected &= has_0p0;
	connected &= has_0m0;
	connected &= has_p00;
	connected &= has_m00;

	return connected;
}

void get_symmops(std::string name_, utility::vector1<core::kinematics::RT> &rt_out, utility::vector1<numeric::xyzVector<core::Real> > &cenop);


int
main( int argc, char * argv [] ) {
try {
	// all spacegroups
	utility::vector1<std::string> spacegroups;


	if (argc<2) {
		spacegroups.push_back("P1");     spacegroups.push_back("P2");      spacegroups.push_back("P21");
		spacegroups.push_back("C2");     spacegroups.push_back("P23");     spacegroups.push_back("F23");
		spacegroups.push_back("I23");    spacegroups.push_back("P213");    spacegroups.push_back("I213");
		spacegroups.push_back("P432");   spacegroups.push_back("P4232");   spacegroups.push_back("F432");
		spacegroups.push_back("F4132");  spacegroups.push_back("I432");    spacegroups.push_back("P4332");
		spacegroups.push_back("P4132");  spacegroups.push_back("I4132");   spacegroups.push_back("P222");
		spacegroups.push_back("P2221");  spacegroups.push_back("P21212");  spacegroups.push_back("P212121");
		spacegroups.push_back("C2221");  spacegroups.push_back("C222");    spacegroups.push_back("F222");
		spacegroups.push_back("I222");   spacegroups.push_back("I212121"); spacegroups.push_back("P4");
		spacegroups.push_back("P41");    spacegroups.push_back("P42");     spacegroups.push_back("P43");
		spacegroups.push_back("I4");     spacegroups.push_back("I41");     spacegroups.push_back("P422");
		spacegroups.push_back("P4212");  spacegroups.push_back("P4122");   spacegroups.push_back("P41212");
		spacegroups.push_back("P4222");  spacegroups.push_back("P42212");  spacegroups.push_back("P4322");
		spacegroups.push_back("P43212"); spacegroups.push_back("I422");    spacegroups.push_back("I4122");
		spacegroups.push_back("P3");     spacegroups.push_back("P31");     spacegroups.push_back("P32");
		spacegroups.push_back("H3");     spacegroups.push_back("P312");    spacegroups.push_back("P321");
		spacegroups.push_back("P3112");  spacegroups.push_back("P3121");   spacegroups.push_back("P3212");
		spacegroups.push_back("P3221");  spacegroups.push_back("H32");     spacegroups.push_back("P6");
		spacegroups.push_back("P61");    spacegroups.push_back("P65");     spacegroups.push_back("P62");
		spacegroups.push_back("P64");    spacegroups.push_back("P63");     spacegroups.push_back("P622");
		spacegroups.push_back("P6122");   spacegroups.push_back("P6522");  spacegroups.push_back("P6222");
		spacegroups.push_back("P6422");   spacegroups.push_back("P6322");
	} else {
		for (int i=1; i<argc; ++i)
			spacegroups.push_back(std::string(argv[i]));
	}

	bool fastmode=false; // run faster for a couple SGs

	for (int sx=1; sx<=(int) spacegroups.size(); ++sx) {
		std::string name=spacegroups[sx];

		std::cerr << "**** " <<  spacegroups[sx] << " ****" << std::endl;

		numeric::xyzMatrix<core::Real> skewM=numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
		if((name=="P3")||(name=="P31")||(name=="P32")||(name=="H3")||(name=="P312")||(name=="P321")||
			(name=="P3112")||(name=="P3121")||(name=="P3212")||(name=="P3221")||(name=="H32")||(name=="P6")||
			(name=="P61")||(name=="P65")||(name=="P62")||(name=="P64")||(name=="P63")||(name=="P622")||
			(name=="P6122")||(name=="P6522")||(name=="P6222")||(name=="P6422")||(name=="P6322"))
			skewM = numeric::xyzMatrix<core::Real>::rows( 1,-0.5,0, 0,0.866025403784439,0, 0,0,1);

		if (name=="P432"||name=="F432"||name=="I432")
			fastmode=true;

		utility::vector1<core::kinematics::RT> rts;
		utility::vector1<numeric::xyzVector<core::Real> > cenops;
		get_symmops( spacegroups[sx], rts, cenops );
		utility::vector1< pointGroupHit > pgs = get_point_groups( rts, skewM, fastmode );

//		for (int i=1; i<=(int)pgs.size(); ++i) {
//			pointGroupHit hit_i = pgs[i];
//			std::cerr << "    " << hit_i.name << ": axis=[" << hit_i.axis[0] <<","<< hit_i.axis[1] <<","<< hit_i.axis[2]
//				<<"]  ori=["<< hit_i.origin[0] <<","<< hit_i.origin[1] <<","<< hit_i.origin[2] <<"] "
//				<<"  lattice_shift=[" << hit_i.ls2[0] <<","<< hit_i.ls2[1] <<","<< hit_i.ls2[2] <<"] "
//				<<" symmops:" << hit_i.symmop2;
//			if (hit_i.symmop3 != 0) std::cerr << " and " << hit_i.symmop3 << std::endl;
//			else std::cerr << std::endl;
//		}

		// sample all pairs of point groups
		utility::vector1< latticeHit > finalHits;
		for (int i=1; i<=(int)pgs.size(); ++i) {
			for (int j=i+1; j<=(int)pgs.size(); ++j) {
				pointGroupHit hit_i = pgs[i];
				pointGroupHit hit_j = pgs[j];

				// add RTs
				utility::vector1<core::kinematics::RT> rts_symm;
				rts_symm.push_back( core::kinematics::RT( rts[hit_i.symmop1].get_rotation(), rts[hit_i.symmop1].get_translation() ) );
				rts_symm.push_back( core::kinematics::RT( rts[hit_i.symmop2].get_rotation(), rts[hit_i.symmop2].get_translation() + hit_i.ls2) );
				if (hit_i.symmop3 > 0)
					rts_symm.push_back( core::kinematics::RT( rts[hit_i.symmop3].get_rotation(), rts[hit_i.symmop3].get_translation() + hit_i.ls3 ) );

				rts_symm.push_back( core::kinematics::RT( rts[hit_j.symmop1].get_rotation(), rts[hit_j.symmop1].get_translation() ) );
				rts_symm.push_back( core::kinematics::RT( rts[hit_j.symmop2].get_rotation(), rts[hit_j.symmop2].get_translation() + hit_j.ls2 ) );
				if (hit_j.symmop3 > 0)   // d symm only
					rts_symm.push_back( core::kinematics::RT( rts[hit_j.symmop3].get_rotation(), rts[hit_j.symmop3].get_translation() + hit_j.ls3 ) );

				bool is_lattice = check_if_forms_lattice( rts_symm, rts, false );   // check if we are fully connected

				if (is_lattice) {
					numeric::xyzVector<core::Real> Iaxis=hit_i.axis, Jaxis=hit_j.axis;

					if (hit_i.name == "D2")
						Iaxis = hit_i.axis.cross( hit_i.axis2 );
					if (hit_j.name == "D2")
						Jaxis = hit_j.axis.cross( hit_j.axis2 );

					Real angle = RAD2DEG * acos( std::max( std::min( Iaxis.dot( Jaxis ), 1.0 ), -1.0) );
					angle = std::min( 180-angle, angle );
					Real offset=0, shift=0;

					if (angle <= 1e-6) {
						offset = ((hit_j.origin-hit_i.origin) - (hit_j.origin-hit_i.origin).dot( Iaxis )* Iaxis).length();
						shift = (hit_j.origin-hit_i.origin).dot( Iaxis );
					} else {
						offset = std::fabs (
							(hit_j.origin-hit_i.origin).dot( hit_i.axis.cross( hit_j.axis ) ) / hit_i.axis.cross( hit_j.axis ).length() );
					}

					bool dup = false;
					for (int q = 1; q<= (int)finalHits.size() && !dup; ++q) {
						if (
							( (finalHits[q].hit1.name == hit_i.name && finalHits[q].hit2.name == hit_j.name) ||
							  (finalHits[q].hit2.name == hit_i.name && finalHits[q].hit1.name == hit_j.name) ) &&
							std::fabs( angle - finalHits[q].angle ) < 1e-6
						) {
							if ((offset==0 && finalHits[q].offset == 0)
								 || (offset!=0 && finalHits[q].offset != 0)) {
								dup = true;
								if (offset<finalHits[q].offset)
									finalHits[q] = latticeHit( hit_i, hit_j, angle, offset, shift );
							}
						}
					}
					if (!dup) {
						finalHits.push_back( latticeHit( hit_i, hit_j, angle, offset, shift ));
					}
				}
			}
		}

		for (int i=1; i<=(int)finalHits.size(); ++i) {
			latticeHit lh = finalHits[i];
			std::cerr << lh.hit1.name << " and " << lh.hit2.name << " at angle = " << lh.angle << " offset = " << lh.offset;
			if (lh.angle == 0) std::cerr << " shift = " << lh.shift;
			std::cerr << std::endl;
			std::cerr <<  "     " <<  lh.hit1.name << " axis=[" << lh.hit1.axis[0] <<","<< lh.hit1.axis[1] <<","<< lh.hit1.axis[2]
				<< "] ";
			if (lh.hit1.name.substr(0,1) == "D")
				std::cerr << " axis2=[" << lh.hit1.axis2[0] <<","<< lh.hit1.axis2[1] <<","<< lh.hit1.axis2[2] << "] ";
			std::cerr << " origin=["<< lh.hit1.origin[0] <<","<< lh.hit1.origin[1] <<","<< lh.hit1.origin[2] <<"]" << std::endl;

			std::cerr <<  "     " <<  lh.hit2.name << " axis=[" << lh.hit2.axis[0] <<","<< lh.hit2.axis[1] <<","<< lh.hit2.axis[2]
				<< "] ";
			if (lh.hit2.name.substr(0,1) == "D")
				std::cerr << " axis2=[" << lh.hit2.axis2[0] <<","<< lh.hit2.axis2[1] <<","<< lh.hit2.axis2[2] << "] ";
			std::cerr << " origin=["<< lh.hit2.origin[0] <<","<< lh.hit2.origin[1] <<","<< lh.hit2.origin[2] <<"]" << std::endl;

		}

		//for (int ss=1; ss<=(int)rts_symm.size(); ++ss)
		//	disp_symmop( rts_symm[ss] );
	}


} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
}
	return 0;
}


///
///
///
void get_symmops(
		std::string name_,
		utility::vector1<core::kinematics::RT> &rt_out,
		utility::vector1<numeric::xyzVector<core::Real> > &cenop) {
	if ( name_ == "P1" ) {
		rt_out.resize(1);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
	}
	if ( name_ == "P121" || name_ == "P2") {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "P1211" || name_ == "P21" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0.5,0) );
	}
	else if ( name_ == "C121" || name_ == "C2" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
	}
	else if ( name_ == "P4" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "P41" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
	}
	else if ( name_ == "P42" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
	}
	else if ( name_ == "P43" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
	}
	else if ( name_ == "I4" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if ( name_ == "I41" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0,0.75) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
 	}
	else if ( name_ == "P422" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "P4212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "P4122" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.25) );
	}
	else if ( name_ == "P41212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
	}
	else if ( name_ == "P4222" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
	}
	else if ( name_ == "P42212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "P4322" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.75) );
	}
	else if ( name_ == "P43212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
	}
	else if ( name_ == "I422" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=8; ++ii) {
			rt_out[8+ii] = rt_out[ii];
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if ( name_ == "I4122" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=8; ++ii) {
			rt_out[8+ii] = rt_out[ii];
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if ( name_ == "P23" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "F23" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.0,0.5,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.0,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.0));
		for (int ii=1; ii<=12; ++ii) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[24+ii] = rt_out[ii];
			rt_out[36+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
	}
	else if ( name_ == "I23" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=12; ++ii) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if ( name_ == "P213" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
	}
	else if ( name_ == "I213" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=12; ++ii) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if ( name_ == "P432" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if ( name_ == "P4232" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
	}
	else if ( name_ == "F432" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.0,0.5,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.0,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.0));
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[48+ii] = rt_out[ii];
			rt_out[72+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
	}
	else if ( name_ == "F4132" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.0,0.5,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.0,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.0));
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[48+ii] = rt_out[ii];
			rt_out[72+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
	}
	else if ( name_ == "I432" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if ( name_ == "P4332" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
	}
	else if ( name_ == "P4132" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.75,0.75,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.75,0.75,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.75) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
	}
	else if ( name_ == "I4132" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if (name_ == "P3") {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P31") {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
	}
	else if (name_ == "P32") {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
	}
	else if (name_ == "H3") {
		rt_out.resize(9);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.666666666666667,0.333333333333333,0.333333333333333));
		cenop.push_back(numeric::xyzVector<Real>(0.333333333333333,0.666666666666667,0.666666666666667));
		for (int ii=1; ii<=3; ++ii) {
			rt_out[3+ii] = rt_out[ii];
			rt_out[3+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
	}
	else if (name_ == "P312") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P321") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P3112") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
	}
	else if (name_ == "P3121") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
	}
	else if (name_ == "P3212") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
	}
	else if (name_ == "P3221") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
	}
	else if (name_ == "H32") {
		rt_out.resize(18);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.666666666666667,0.333333333333333,0.333333333333333));
		cenop.push_back(numeric::xyzVector<Real>(0.333333333333333,0.666666666666667,0.666666666666667));
		for (int ii=1; ii<=6; ++ii) {
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
	}
	else if (name_ == "P6") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P61") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
	}
	else if (name_ == "P65") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
	}
	else if (name_ == "P62") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
	}
	else if (name_ == "P64") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
	}
	else if (name_ == "P63") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
	}
	else if (name_ == "P622") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P6122") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
	}
	else if (name_ == "P6522") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
	}
	else if (name_ == "P6222") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
	}
	else if (name_ == "P6422") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
	}
	else if (name_ == "P6322") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P222") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
	}
	else if (name_ == "P2221") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
	}
	else if (name_ == "P21212") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
	}
	else if (name_ == "P212121") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
	}
	else if (name_ == "C2221") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
	}
	else if (name_ == "C222") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
	}
	else if (name_ == "F222") {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0,0.5,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0,0.5));
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[8+ii] = rt_out[ii];
			rt_out[12+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
	}
	else if (name_ == "I222") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
	else if (name_ == "I212121") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		// cenops
		cenop.push_back(numeric::xyzVector<Real>(0.5,0.5,0.5));
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
	}
}
