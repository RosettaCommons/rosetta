// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file
/// @brief

#include <apps/pilot/frank/spacegroup.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RB_geometry.hh>

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

static THREAD_LOCAL basic::Tracer TR( "cryst.design" );

OPT_1GRP_KEY(String, crystdock, spacegroup)
OPT_1GRP_KEY(Real, crystdock, A)
OPT_1GRP_KEY(Real, crystdock, B)
OPT_1GRP_KEY(Real, crystdock, C)
OPT_1GRP_KEY(Real, crystdock, alpha)
OPT_1GRP_KEY(Real, crystdock, beta)
OPT_1GRP_KEY(Real, crystdock, gamma)
OPT_1GRP_KEY(Real, crystdock, maxclash)
OPT_1GRP_KEY(Real, crystdock, mininterface)
OPT_1GRP_KEY(Real, crystdock, mininterfacesum)
OPT_1GRP_KEY(Real, crystdock, trans_step)
OPT_1GRP_KEY(Real, crystdock, rot_step)
OPT_1GRP_KEY(Real, crystdock, cb_radius)
OPT_1GRP_KEY(Integer, crystdock, nmodels)
OPT_1GRP_KEY(Integer, crystdock, rotnum)
OPT_1GRP_KEY(Integer, crystdock, symnum)
OPT_1GRP_KEY(Boolean, crystdock, ssonly)
OPT_1GRP_KEY(Boolean, crystdock, debug)
OPT_1GRP_KEY(Boolean, crystdock, debug_exact)
OPT_1GRP_KEY(Boolean, crystdock, eval_native)
OPT_1GRP_KEY(Boolean, crystdock, randomize_orientation)
OPT_1GRP_KEY(Real, crystdock, n_clashdist)
OPT_1GRP_KEY(Real, crystdock, ca_clashdist)
OPT_1GRP_KEY(Real, crystdock, c_clashdist)
OPT_1GRP_KEY(Real, crystdock, o_clashdist)
OPT_1GRP_KEY(Real, crystdock, cb_clashdist)
OPT_1GRP_KEY(Real, crystdock, sigwidth)
OPT_1GRP_KEY(Real, crystdock, interfacedist)
OPT_1GRP_KEY(Real, crystdock, interface_sigwidth)
OPT_1GRP_KEY(Real, crystdock, cluster_cutoff)
OPT_1GRP_KEY(Integer, crystdock, expand_rounds)
OPT_1GRP_KEY(Real, crystdock, random_rotate)
OPT_1GRP_KEY(Boolean, crystdock, compact)
OPT_KEY(Integer,  run_i)
OPT_KEY(Integer,  run_j)

////////////////////////////////////////////////
// helper functions
inline int pos_mod(int x,int y) {
	int r=x%y; if ( r<0 ) r+=y;
	return r;
}
inline Real pos_mod(Real x,Real y) {
	Real r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
inline int min_mod(int x,int y) {
	int r=x%y; if ( r<-y/2 ) { r+=y; } if ( r>=y/2 ) { r-=y; }
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y); if ( r<-0.5*y ) { r+=y; } if ( r>=0.5*y ) { r-=y; }
	return r;
}

inline bool transforms_equiv(
	numeric::xyzMatrix<Real> const &S1, numeric::xyzVector<Real> const &T1,
	numeric::xyzMatrix<Real> const &S2, numeric::xyzVector<Real> const &T2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() ) +
		std::fabs( T1[0] - T2[0] ) +  std::fabs( T1[1] - T2[1] ) +  std::fabs( T1[2] - T2[2] );
	return (err <= 1e-6);
}

inline double sign(double x) {return (x >= 0.0) ? 1.0 : -1.0;}
inline double norm4(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

/// backward interpolation (unrotated->rotated)
///    used to generate grids
template <class S>
core::Real interp_linear(
	ObjexxFCL::FArray3D< S > const & data ,
	numeric::xyzVector< core::Real > const & idxX ) {

	numeric::xyzVector<int> pt000, pt111;
	numeric::xyzVector< core::Real > fpart,neg_fpart;
	numeric::xyzVector<int> srcgrid(data.u1(),data.u2(),data.u3());

	pt000[0] = (int)(std::floor(idxX[0]+1e-6));
	pt000[1] = (int)(std::floor(idxX[1]+1e-6));
	pt000[2] = (int)(std::floor(idxX[2]+1e-6));

	// interpolation coeffs
	fpart[0] = idxX[0]-pt000[0]; neg_fpart[0] = 1-fpart[0];
	fpart[1] = idxX[1]-pt000[1]; neg_fpart[1] = 1-fpart[1];
	fpart[2] = idxX[2]-pt000[2]; neg_fpart[2] = 1-fpart[2];
	S retval = (S)0.0;

	// bound check
	if ( pt000[0] < -srcgrid[0]/2 || pt000[0] > srcgrid[0]/2 ) return 0.0;
	if ( pt000[1] < -srcgrid[1]/2 || pt000[1] > srcgrid[1]/2 ) return 0.0;
	if ( pt000[2] < -srcgrid[2]/2 || pt000[2] > srcgrid[2]/2 ) return 0.0;
	if ( pt000[0] < 0 ) pt000[0] += srcgrid[0];
	if ( pt000[1] < 0 ) pt000[1] += srcgrid[1];
	if ( pt000[2] < 0 ) pt000[2] += srcgrid[2];

	pt111[0] = (pt000[0]+1); if ( pt111[0]>=srcgrid[0] ) pt111[0]=0;
	pt111[1] = (pt000[1]+1); if ( pt111[1]>=srcgrid[1] ) pt111[1]=0;
	pt111[2] = (pt000[2]+1); if ( pt111[2]>=srcgrid[2] ) pt111[2]=0;

	retval += neg_fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt000[0]+1,pt000[1]+1,pt000[2]+1);
	retval += neg_fpart[0]*neg_fpart[1]*    fpart[2] * data(pt000[0]+1,pt000[1]+1,pt111[2]+1);
	retval += neg_fpart[0]*    fpart[1]*neg_fpart[2] * data(pt000[0]+1,pt111[1]+1,pt000[2]+1);
	retval += neg_fpart[0]*    fpart[1]*    fpart[2] * data(pt000[0]+1,pt111[1]+1,pt111[2]+1);
	retval +=     fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt111[0]+1,pt000[1]+1,pt000[2]+1);
	retval +=     fpart[0]*neg_fpart[1]*    fpart[2] * data(pt111[0]+1,pt000[1]+1,pt111[2]+1);
	retval +=     fpart[0]*    fpart[1]*neg_fpart[2] * data(pt111[0]+1,pt111[1]+1,pt000[2]+1);
	retval +=     fpart[0]*    fpart[1]*    fpart[2] * data(pt111[0]+1,pt111[1]+1,pt111[2]+1);
	return (core::Real)retval;
}


// rotation stuff
void
euler2rot( core::Real a, core::Real b, core::Real g, numeric::xyzMatrix<Real> &R) {
	R.xx() = -sin(DEG2RAD*a)*cos(DEG2RAD*b)*sin(DEG2RAD*g) + cos(DEG2RAD*a)*cos(DEG2RAD*g);
	R.xy() =  cos(DEG2RAD*a)*cos(DEG2RAD*b)*sin(DEG2RAD*g) + sin(DEG2RAD*a)*cos(DEG2RAD*g);
	R.xz() =  sin(DEG2RAD*b)*sin(DEG2RAD*g);
	R.yx() = -sin(DEG2RAD*a)*cos(DEG2RAD*b)*cos(DEG2RAD*g) - cos(DEG2RAD*a)*sin(DEG2RAD*g);
	R.yy() =  cos(DEG2RAD*a)*cos(DEG2RAD*b)*cos(DEG2RAD*g) - sin(DEG2RAD*a)*sin(DEG2RAD*g);
	R.yz() =  sin(DEG2RAD*b)*cos(DEG2RAD*g);
	R.zx() =  sin(DEG2RAD*a)*sin(DEG2RAD*b);
	R.zy() = -cos(DEG2RAD*a)*sin(DEG2RAD*b);
	R.zz() =  cos(DEG2RAD*b);
}

// angle of rotation from rotation matrix
// an abomination of 0- and 1- indexing
core::Real
R2ang(numeric::xyzMatrix<Real> R) {
	Real q0 = ( R(1,1) + R(2,2) + R(3,3) + 1.0) / 4.0;
	Real q1 = ( R(1,1) - R(2,2) - R(3,3) + 1.0) / 4.0;
	Real q2 = (-R(1,1) + R(2,2) - R(3,3) + 1.0) / 4.0;
	Real q3 = (-R(1,1) - R(2,2) + R(3,3) + 1.0) / 4.0;

	if ( q0 < 0.0 ) q0 = 0.0f; else q0 = sqrt(q0);
	if ( q1 < 0.0 ) q1 = 0.0f; else q1 = sqrt(q1);
	if ( q2 < 0.0 ) q2 = 0.0f; else q2 = sqrt(q2);
	if ( q3 < 0.0 ) q3 = 0.0f; else q3 = sqrt(q3);

	if ( q0 >= q1 && q0 >= q2 && q0 >= q3 ) {
		q1 *= sign(R(3,2) - R(2,3));
		q2 *= sign(R(1,3) - R(3,1));
		q3 *= sign(R(2,1) - R(1,2));
	} else if ( q1 >= q0 && q1 >= q2 && q1 >= q3 ) {
		q0 *= sign(R(3,2) - R(2,3));
		q2 *= sign(R(2,1) + R(1,2));
		q3 *= sign(R(1,3) + R(3,1));
	} else if ( q2 >= q0 && q2 >= q1 && q2 >= q3 ) {
		q0 *= sign(R(1,3) - R(3,1));
		q1 *= sign(R(2,1) + R(1,2));
		q3 *= sign(R(3,2) + R(2,3));
	} else if ( q3 >= q0 && q3 >= q1 && q3 >= q2 ) {
		q0 *= sign(R(2,1) - R(1,2));
		q1 *= sign(R(3,1) + R(1,3));
		q2 *= sign(R(3,2) + R(2,3));
	}
	Real r = norm4(q0, q1, q2, q3);
	q0 /= r; q1 /= r; q2 /= r; q3 /= r;
	Real angle = 2*RAD2DEG*std::fabs(acos(q0));

	return (angle);
}


////////////////////////////////////////////////
// Uniformly sample rotation space
//   * use idea from A. Yershova, et al, International Journal of Robotics Research, 2009.
//       to get SO3 samples from S1+S2 uniform sampling
//   * use icosahedral embedding to generate approximately uniform coverage of S2
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
		if ( R.xx() + R.yy() + R.zz() > 0 ) {
			S = 0.5 / sqrt(R.xx() + R.yy() + R.zz());
			w_ = 0.25 / S;
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
		if ( 1-w_*w_ < 1e-6 ) {
			return numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);
		}

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


class UniformRotationSampler {
private:
	utility::vector1< Quat > rotlist_;
	Real theta_;

public:
	Size
	nrots() const { return rotlist_.size(); }

	void get(Size ii, numeric::xyzMatrix<Real> &R) const {
		if ( ii==0 ) {
			R = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
		} else {
			R = rotlist_[ii].asR();
		}
	}

	void
	generateIcosahedralSamples(utility::vector1<numeric::xyzVector<Real> > &ico, Size nsub) {
		using numeric::constants::d::pi;
		Real theta = 26.56505117707799 * pi / 180.0;
		Real stheta = sin(theta);
		Real ctheta = cos(theta);
		ico.clear();

		ico.push_back(numeric::xyzVector<Real>(0.0,0.0,-1.0));
		Real phi = pi / 5.0;
		for ( int i=0; i<5; ++i ) {
			ico.push_back(numeric::xyzVector<Real>( ctheta * cos(phi), ctheta * sin(phi), -stheta ));
			phi += 2.0 * pi / 5.0;
		}
		phi = 0.0;
		for ( int i=0; i<5; ++i ) {
			ico.push_back(numeric::xyzVector<Real>( ctheta * cos(phi), ctheta * sin(phi), stheta ));
			phi += 2.0 * pi / 5.0;
		}
		ico.push_back(numeric::xyzVector<Real>(0.0,0.0,1.0));
		Size TRIS[][3] = {
			{0,2,1}, {0,3,2}, {0,4,3}, {0,5,4}, {0,1,5}, {1,2,7}, {2,3,8}, {3,4,9}, {4,5,10}, {5,1,6},
			{1,7,6}, {2,8,7}, {3,9,8}, {4,10,9}, {5,6,10}, {6,7,11}, {7,8,11}, {8,9,11}, {9,10,11}, {10,6,11}
			};
		Size EDGES[][2] =  {
			{0,1}, {0,2}, {0,3}, {0,4}, {0,5}, {1,2}, {1,5}, {1,6}, {1,7}, {2,3}, {2,7}, {2,8}, {3,4}, {3,8}, {3,9},
			{4,5}, {4,9}, {4,10}, {5,6}, {5,10}, {6,7}, {6,10}, {6,11}, {7,8}, {7,11}, {8,9}, {8,11}, {9,10}, {9,11}, {10,11}
			};

		// subdivide edges
		for ( int i=0; i<30; ++i ) {
			numeric::xyzVector<Real> a = ico[EDGES[i][0]+1];
			numeric::xyzVector<Real> b = ico[EDGES[i][1]+1];
			numeric::xyzVector<Real> ab = b-a;
			for ( int j=1; j<(int)nsub; ++j ) {
				numeric::xyzVector<Real> new_j;
				new_j = a+(((Real)j)/((Real)nsub))*ab;
				new_j = new_j/new_j.length();
				ico.push_back(new_j);
			}
		}

		// subdivide faces
		for ( int i=0; i<20; ++i ) {
			numeric::xyzVector<Real> a = ico[TRIS[i][0]+1];
			numeric::xyzVector<Real> b = ico[TRIS[i][1]+1];
			numeric::xyzVector<Real> c = ico[TRIS[i][2]+1];
			numeric::xyzVector<Real> ab = b-a;
			numeric::xyzVector<Real> ac = c-a;
			for ( int j=1; j<(int)nsub-1; ++j ) {
				for ( int k=j; k<(int)nsub-1; ++k ) {
					numeric::xyzVector<Real> new_j;
					new_j = a + ((Real)j)/((Real)nsub)*ac + ((Real)(nsub-1-k))/((Real)nsub)*ab;
					new_j = new_j/new_j.length();
					ico.push_back(new_j);
				}
			}
		}
	}

	UniformRotationSampler(Real theta) {
		theta_ = theta;
		if ( theta<=1e-6 ) {
			rotlist_.resize(1);
			return;
		}

		using numeric::constants::d::pi;
		int S1subs = (int)std::floor(360.0 / theta + 0.5);
		int S2subs = (int)std::floor(60.0 / theta + 0.5);

		utility::vector1<numeric::xyzVector<Real> > ico;
		generateIcosahedralSamples(ico, S2subs);

		int nsamples = S1subs*ico.size();
		rotlist_.resize( nsamples );
		int counter=1;

		// precompute cos/sin (psi/2)for S1
		utility::vector1<core::Real> cos_psi_over_2(S1subs), sin_psi_over_2(S1subs);
		for ( int j=0; j<S1subs; ++j ) {
			Real psi_over_2 = j*pi/((Real)S1subs);
			cos_psi_over_2[j+1] = cos(psi_over_2);
			sin_psi_over_2[j+1] = sin(psi_over_2);
		}
		for ( int i=1; i<=(int)ico.size(); ++i ) {
			// convert to spherical ASSUMES NORMALIZED
			Real theta = pi-acos( ico[i][2] );
			Real phi = atan2( ico[i][1], ico[i][0] ) + pi;
			for ( int j=0; j<S1subs; ++j ) {
				Real psi = 2.0*j*pi/((Real)S1subs);
				Real w = sin(theta/2)*sin(phi+psi/2);
				Real x = cos(theta/2)*cos(psi/2);
				Real y = cos(theta/2)*sin(psi/2);
				Real z = sin(theta/2)*cos(phi+psi/2);
				rotlist_[counter++] = Quat(x,y,z,w);
			}
		}
	}

	void
	remove_redundant(
		utility::vector1<core::kinematics::RT> const &symmops,
		numeric::xyzMatrix<Real> const &i2c, numeric::xyzMatrix<Real> const &c2i ) {

		utility::vector1<bool> tokeep(rotlist_.size(), true);
		utility::vector1<numeric::xyzMatrix<Real> > Rs(rotlist_.size());
		utility::vector1<numeric::xyzMatrix<Real> > Rinvs(rotlist_.size());

		for ( int i=1; i<=(int)rotlist_.size(); ++i ) {
			Rs[i] = rotlist_[i].asR();
			Rinvs[i] = numeric::inverse( Rs[i] );
		}


		for ( int s=2; s<=(int)symmops.size(); ++s ) {
			numeric::xyzMatrix<Real> S = symmops[s].get_rotation();
			if ( S.xy()==0 && S.xz()==0 && S.yz()==0 && S.yx()==0 && S.zx()==0 && S.zy()==0 && S.xx()==1 && S.yy()==1 && S.zz()==1 ) {
				continue; // identity
			}

			// cartesian-space S
			numeric::xyzMatrix<Real> Scart = i2c*S*c2i;

			for ( int i=1; i<=(int)rotlist_.size(); ++i ) {
				if ( !tokeep[i] ) { continue; }
				numeric::xyzMatrix<Real> SR = Scart*Rs[i];
				for ( int j=i+1; j<=(int)rotlist_.size(); ++j ) {
					if ( !tokeep[j] ) { continue; }

					Real ang_ij = R2ang( SR*Rinvs[j] );
					if ( ang_ij<0.65*theta_ ) {
						tokeep[j]=false;
					}
				}
			}
		}

		utility::vector1< Quat > rotlist_old = rotlist_;
		rotlist_.clear();
		for ( int i=1; i<=(int)rotlist_old.size(); ++i ) {
			if ( tokeep[i] ) rotlist_.push_back(rotlist_old[i]);
		}

		TR << "Trimmed rot list has " << rotlist_.size() << " rotations (was " << rotlist_old.size() << ")\n";
	}
};


///////////////////////////////////////////////////////////////////////////////
struct SingleInterface {
	SingleInterface( numeric::xyzMatrix<core::Real> R_in, numeric::xyzVector<core::Real> T_in, Real cb_overlap_in ) {
		cb_overlap_ = cb_overlap_in;
		R_ = R_in;
		T_ = T_in;
	}

	Real cb_overlap_;
	numeric::xyzMatrix<core::Real> R_;
	numeric::xyzVector<core::Real> T_;
};

///////////////////////////////////////////////////////////////////////////////

// fpd this class is used for storing hits over all rotations
//    same top-N priority_queue but much larger store
struct InterfaceHit {
	InterfaceHit( Real score_in, Real x_in, Real y_in, Real z_in, Size rot_index_in ) {
		score=score_in;
		x=x_in; y=y_in; z=z_in;
		rot_index=rot_index_in;
	}
	Real score;
	Real x,y,z;
	Size rot_index;
	//utility::vector1<SingleInterface> iinfo;

	utility::vector1< std::string >
	to_string() {
		utility::vector1< std::string > retval;
		std::ostringstream oss;
		oss << "[rot " << rot_index << "] score = "<< score;
		retval.push_back( oss.str() );

		return retval;
	}
};

class InterfaceHitComparitor {
public:
	bool operator()(InterfaceHit& t1, InterfaceHit& t2) {
		return (t1.score > t2.score);
	}
};

class InterfaceHitDatabase  {
private:
	core::Size N_;
	std::priority_queue<InterfaceHit, std::vector<InterfaceHit>, InterfaceHitComparitor> queue_;

public:
	InterfaceHitDatabase(int N_in) { N_ = N_in; }

	void add_interface( InterfaceHit interface ) {
		if ( queue_.size() <N_ ) {
			queue_.push( interface );
		} else if ( interface.score > queue_.top().score ) {
			queue_.pop();
			queue_.push( interface );
		}
	}

	InterfaceHit pop() {
		InterfaceHit retval = queue_.top();
		queue_.pop();
		return retval;
	}

	InterfaceHit top() {
		InterfaceHit retval = queue_.top();
		return retval;
	}

	Size size() { return queue_.size(); }
};


///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////

class CrystDock : public protocols::moves::Mover {
private:
	core::Real ca_clashdist_, cb_clashdist_, n_clashdist_, c_clashdist_, o_clashdist_;
	core::Real interfacedist_,  voxel_volume_;
	numeric::xyzVector<int> grid_, oversamplegrid_;
	numeric::xyzMatrix<Real> i2c_, c2i_;
	Spacegroup sg_;

	// parameters from options
	core::Real maxclash_, mininterface_, trans_step_, rot_step_;
	Size nmodels_, rotnum_;
	bool ss_only_, eval_native_;
	bool debug_, debug_exact_;
	bool compact_;

public:
	CrystDock() {
		// to do: make these options
		n_clashdist_  =  option[crystdock::n_clashdist];   // ros VDW=1.75
		ca_clashdist_ =  option[crystdock::ca_clashdist];   // ros VDW=2.0
		c_clashdist_  =  option[crystdock::c_clashdist];   // ros VDW=2.0
		o_clashdist_  =  option[crystdock::o_clashdist];   // ros VDW=1.55
		cb_clashdist_ =  option[crystdock::cb_clashdist];   // ros VDW=2.0    .. make this a bit smaller to encourage better packing
		interfacedist_ = option[crystdock::interfacedist];

		mininterface_ = option[ crystdock::mininterface ]();
		maxclash_ = option[ crystdock::maxclash ]();
		trans_step_ = option[ crystdock::trans_step ]();
		rot_step_ = option[ crystdock::rot_step ]();
		nmodels_ = option[ crystdock::nmodels ]();
		rotnum_ = option[ crystdock::rotnum ]();
		debug_ = option[ crystdock::debug ]();
		debug_exact_ = option[ crystdock::debug_exact ]();
		ss_only_ = option[ crystdock::ssonly ]();
		eval_native_ = option[ crystdock::eval_native ]();
		compact_=option[ crystdock::compact ]();
	}

	virtual std::string get_name() const { return "CrystDock"; }

	// write density grids in MRC format for debugging
	void
	writeMRC(FArray3D<Real> density, std::string mapfilename, bool is_oversampled=false, bool fftshift=false);

	// build occupancy and interaction masks from pose
	// sample in bounding box around molecule (rather than unit cell)
	// TODO: first align pose on principal axes
	void setup_maps( Pose & pose, FArray3D<Real> &rho_ca, FArray3D<Real> &rho_cb, Real trans_step);

	// resample maps subject to rotation
	// get self clashes and self rotations
	core::Real resample_maps_and_get_self(
		FArray3D<Real> const &rho_ca, FArray3D<Real> const &rho_cb,
		numeric::xyzMatrix<Real> R, Spacegroup const &sg,
		FArray3D<Real> &r_rho_ca, FArray3D<Real> &r_rho_cb,
		utility::vector1<SingleInterface> &p1_interface_map );

	// recenter pose
	numeric::xyzVector<Real> center_pose_at_origin( Pose & pose );

	// add CRYST1 line
	void add_crystinfo_to_pose( Pose & pose );

	// aplpy xform and dump PDB
	void dump_transformed_pdb( Pose pose, InterfaceHit ih, UniformRotationSampler const &urs, std::string outname,std::string basename  );

	// nearest-neighbor interpolation subject to grid-space transform
	void transform_map(
		FArray3D<Real> const &rho,
		numeric::xyzMatrix<Real> S, numeric::xyzVector<Real> T,
		FArray3D<Real> &Srho);

	// same as previous function, assumes offset of 0
	void transform_map_offset0(
		FArray3D<Real> const &rho,
		numeric::xyzMatrix<Real> S,
		FArray3D<Real> &Srho);

	// get maximum radius from a pose
	Real get_radius_of_pose( Pose & pose );

	//Calculate transformed distance
	Real get_transform_distance (InterfaceHit ih_vec, InterfaceHit ih_vec_clustered , UniformRotationSampler const &urs, Real radius);
	// wrapper does a series of 2d ffts
	//    this could be an option to the low-level fft code to possibly save time
	void
	fft2dslice( FArray3D<Real> const &rho, FArray3D< std::complex<Real> > &Frho, int axisSlice );

	// wrapper does a series of 2d iffts
	void
	ifft2dslice( FArray3D< std::complex<Real> > const &Frho, FArray3D<Real> &rho, int axisSlice );

	// project scores along some axis, not necessarily orthogonal to unit cell
	void
	project_along_axis( FArray3D<Real> &rho, numeric::xyzVector<Real> axis );

	// the main convolution operation used for both CA and CB maps
	void
	do_convolution( FArray3D<Real> const &rho, FArray3D<Real> const &Srho, numeric::xyzMatrix<Real> S, FArray3D<Real> &conv_out);

	// get the score per interface
	core::Real
	get_interface_score(
		Pose pose,
		utility::vector1<SingleInterface> &allInterfaces,
		utility::vector1<core::kinematics::RT> const &rts_all );

	// get the score per interface
	Real
	get_interfaces(
		utility::vector1<core::kinematics::RT> const &rts,
		numeric::xyzMatrix<Real> R,
		numeric::xyzVector<Real> xyz_grid,
		utility::vector1<Size> symmopList,
		FArray3D<Real> const &rho_cb,
		utility::vector1<SingleInterface> &allInterfaces );

	// get the score per interface
	void
	get_interfaces_allatom(
		Pose pose,
		utility::vector1<core::kinematics::RT> const &rts,
		numeric::xyzMatrix<Real> R,
		numeric::xyzVector<Real> xyz_grid,
		utility::vector1<SingleInterface> &allInterfaces );

	// get the exact clash score (debugging/eval_native only)
	core::Real
	get_clash_score_exact(
		numeric::xyzVector<int> xyz_grid,
		numeric::xyzMatrix<Real> R,
		numeric::xyzVector<Real> T,
		FArray3D<Real> const &r_rho_ca);

	void apply( Pose & pose);
};

// write density grids in MRC --  debugging only for now
void
CrystDock::writeMRC(FArray3D<Real> density, std::string mapfilename, bool is_oversampled /*=false*/, bool fftshift /*=false*/) {
	const int CCP4HDSIZE = 1024;  // size of CCP4/MRC header
	std::fstream outx( (mapfilename).c_str() , std::ios::binary | std::ios::out );

	float buff_f, buff_vf[3];
	int buff_i, buff_vi[3], symBytes = 0;

	if ( !outx ) {
		TR << "Error opening " << mapfilename << " for writing." << std::endl;
		return;
	}

	// extent
	buff_vi[0] = density.u1(); buff_vi[1] = density.u2(); buff_vi[2] = density.u3();
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// mode
	buff_i = 2;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// origin
	int ori_int[3] = {0,0,0};
	outx.write(reinterpret_cast <char*>(ori_int), sizeof(int)*3);

	// grid
	int grid[3] = {density.u1(),density.u2(),density.u3()};
	outx.write(reinterpret_cast <char*>(&grid[0]), sizeof(int)*3);

	// cell params
	Real Aeff=sg_.A(), Beff=sg_.B(), Ceff=sg_.C();
	if ( is_oversampled ) {
		Aeff *= ((Real)oversamplegrid_[0]) / ((Real)grid_[0]);
		Beff *= ((Real)oversamplegrid_[1]) / ((Real)grid_[1]);
		Ceff *= ((Real)oversamplegrid_[2]) / ((Real)grid_[2]);
	}

	float cellDimensions[3] = {(float)Aeff,(float)Beff,(float)Ceff};
	float cellAngles[3] = {(float)sg_.alpha(),(float)sg_.beta(),(float)sg_.gamma()};
	outx.write(reinterpret_cast <char*>(&cellDimensions), sizeof(float)*3);
	outx.write(reinterpret_cast <char*>(&cellAngles), sizeof(float)*3);

	// crs2xyz
	buff_vi[0] = 1; buff_vi[1] = 2; buff_vi[2] = 3;
	outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

	// min, max, mean dens
	buff_vf[0] = -100.0; buff_vf[1] = 100.0; buff_vf[2] = 0.0;
	outx.write(reinterpret_cast <char*>(buff_vf), sizeof(float)*3);

	// 4 bytes junk
	buff_i = 0;
	outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

	// symmops
	outx.write(reinterpret_cast <char*>(&symBytes), sizeof(int));

	// 104 bytes junk
	buff_i = 0;
	for ( int i=0; i<25; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// alt origin (MRC)
	float ori_float[3]={0,0,0};
	outx.write(reinterpret_cast <char*>(ori_float), sizeof(float)*3);

	// Write "MAP" at byte 208, indicating a CCP4 file.
	char buff_s[80]; strcpy(buff_s, "MAP DD");
	outx.write(reinterpret_cast <char*>(buff_s), 8);

	// fill remainder of head with junk
	int nJunkWords = (CCP4HDSIZE - 216) /4;
	buff_i = 0;
	for ( int i=0; i<nJunkWords; ++i ) {
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
	}

	// data
	if ( fftshift ) {
		int coord[3], rcoord[3];
		for ( coord[2] = 1; coord[2] <= density.u3(); coord[2]++ ) {
			for ( coord[1] = 1; coord[1] <= density.u2(); coord[1]++ ) {
				for ( coord[0] = 1; coord[0] <= density.u1(); coord[0]++ ) {
					rcoord[0] = pos_mod( coord[0]-1+density.u1()/2, density.u1())+1;
					rcoord[1] = pos_mod( coord[1]-1+density.u2()/2, density.u2())+1;
					rcoord[2] = pos_mod( coord[2]-1+density.u3()/2, density.u3())+1;
					buff_f = (float) density(rcoord[0],rcoord[1],rcoord[2]);
					outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
				}
			}
		}
	} else {
		int coord[3];
		for ( coord[2] = 1; coord[2] <= density.u3(); coord[2]++ ) {
			for ( coord[1] = 1; coord[1] <= density.u2(); coord[1]++ ) {
				for ( coord[0] = 1; coord[0] <= density.u1(); coord[0]++ ) {
					buff_f = (float) density(coord[0],coord[1],coord[2]);
					outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
				}
			}
		}
	}
}


// build occupancy and interaction masks from pose
// since we will be sampling rotations from this map, we sample over a larger volume than the unit cell
//  ASSUMES POSE IS CENTERED
void
CrystDock::setup_maps( Pose & pose, FArray3D<Real> &rho_ca, FArray3D<Real> &rho_cb, Real trans_step) {
	Real ATOM_MASK_PADDING = 2.0;
	Real UNIT_CELL_PADDING = 6.0;  // IN GRID POINTS!
	Real sigwidth=option[crystdock::sigwidth];
	Real interface_sigwidth=option[crystdock::interface_sigwidth];
	Size minmult = sg_.minmult();

	// find true grid
	grid_[0] = minmult*(Size)std::floor(sg_.A()/(minmult*trans_step) + 0.5);
	grid_[1] = minmult*(Size)std::floor(sg_.B()/(minmult*trans_step) + 0.5);
	grid_[2] = minmult*(Size)std::floor(sg_.C()/(minmult*trans_step) + 0.5);

	numeric::xyzMatrix<Real> i2f, f2i;
	f2i = numeric::xyzMatrix<Real>::rows( (Real)grid_[0],0,0,  0,(Real)grid_[1],0,  0,0,(Real)grid_[2] );
	i2f = numeric::xyzMatrix<Real>::rows( 1.0/((Real)grid_[0]),0,0,  0,1.0/((Real)grid_[1]),0,  0,0,1.0/((Real)grid_[2]) );

	i2c_ = sg_.f2c()*i2f;
	c2i_ = f2i*sg_.c2f();

	// find oversampled grid
	oversamplegrid_ = numeric::xyzVector<int>(0,0,0);
	for ( int i=1 ; i<=(int)pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		for ( int j=1; j<=4; ++j ) {
			numeric::xyzVector< core::Real> xyz_j = pose.residue(i).atom(j).xyz();
			numeric::xyzVector< core::Real> atm_idx = c2i_*xyz_j;
			for ( int k=0; k<3; ++k ) {
				oversamplegrid_[k] = std::max(oversamplegrid_[k], 2*(int)std::floor( (std::abs(atm_idx[k])+UNIT_CELL_PADDING) )+1 );
			}
		}
	}
	TR << "Base grid = [" <<  oversamplegrid_[0] << " , " <<  oversamplegrid_[1]<< " , " << oversamplegrid_[2] << "]" << std::endl;

	rho_ca.dimension( oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2] ); rho_ca = 1;
	rho_cb.dimension( oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2] ); rho_cb = 1;

	// loop over bb heavyatoms
	for ( int i=1 ; i<=(int)pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;

		for ( int j=1; j<=4; ++j ) {
			core::Real clashdist=0.0;
			if ( j==1 ) clashdist = n_clashdist_;
			if ( j==2 ) clashdist = ca_clashdist_;
			if ( j==3 ) clashdist = c_clashdist_;
			if ( j==4 ) clashdist = o_clashdist_;

			numeric::xyzVector< core::Real> xyz_j = pose.residue(i).atom(j).xyz();
			numeric::xyzVector< core::Real> atm_idx = c2i_*xyz_j;
			numeric::xyzVector< core::Real> atm_j, del_ij;

			for ( int z=1; z<=oversamplegrid_[2]; ++z ) {
				atm_j[2] = z-1;
				del_ij[2] = min_mod(atm_idx[2] - atm_j[2], (Real)oversamplegrid_[2]);
				del_ij[0] = del_ij[1] = 0.0;
				if ( (i2c_*del_ij).length_squared() > (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING) ) continue;
				for ( int y=1; y<=oversamplegrid_[1]; ++y ) {
					atm_j[1] = y-1;
					del_ij[1] = min_mod(atm_idx[1] - atm_j[1], (Real)oversamplegrid_[1]);
					del_ij[0] = 0.0;
					if ( (i2c_*del_ij).length_squared() > (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING) ) continue;
					for ( int x=1; x<=oversamplegrid_[0]; ++x ) {
						atm_j[0] = x-1;
						del_ij[0] = min_mod(atm_idx[0] - atm_j[0], (Real)oversamplegrid_[0]);
						numeric::xyzVector< core::Real > cart_del_ij = i2c_*del_ij;
						core::Real d2 = (cart_del_ij).length_squared();
						if ( d2 <= (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING) ) {
							core::Real doff = sqrt(d2) - clashdist;
							core::Real sig = 1 / ( 1 + exp ( -sigwidth*doff ) );   //
							rho_ca(x,y,z) *= sig;
						}
					}
				}
			}
		}
	}

	// loop over CBs
	for ( int i=1 ; i<=(int)pose.total_residue(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_gly ) continue;
		numeric::xyzVector< core::Real> CB = pose.residue(i).atom(5).xyz();
		numeric::xyzVector< core::Real> atm_idx = c2i_*CB;
		numeric::xyzVector< core::Real> atm_j, del_ij;

		for ( int z=1; z<=oversamplegrid_[2]; ++z ) {
			atm_j[2] = z-1;
			del_ij[2] = min_mod(atm_idx[2] - atm_j[2], (Real)oversamplegrid_[2]);
			del_ij[0] = del_ij[1] = 0.0;
			if ( (i2c_*del_ij).length_squared() > (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING) ) continue;
			for ( int y=1; y<=oversamplegrid_[1]; ++y ) {
				atm_j[1] = y-1;
				del_ij[1] = min_mod(atm_idx[1] - atm_j[1], (Real)oversamplegrid_[1]);
				del_ij[0] = 0.0;
				if ( (i2c_*del_ij).length_squared() > (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING) ) continue;
				for ( int x=1; x<=oversamplegrid_[0]; ++x ) {
					atm_j[0] = x-1;
					del_ij[0] = min_mod(atm_idx[0] - atm_j[0], (Real)oversamplegrid_[0]);
					numeric::xyzVector< core::Real > cart_del_ij = i2c_*del_ij;
					core::Real d2 = (cart_del_ij).length_squared();
					if ( d2 <= (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING) ) {
						core::Real doff = sqrt(d2) - cb_clashdist_;
						core::Real sig = 1 / ( 1 + exp ( -sigwidth*doff ) );   // '6' gives sigmoid dropoff
						rho_ca(x,y,z) *= sig;

						if ( !ss_only_ || pose.secstruct(i)!='L' ) {
							doff = sqrt(d2) - interfacedist_;
							sig = 1 / ( 1 + exp ( -interface_sigwidth*doff ) );
							rho_cb(x,y,z) *= sig;
						}
					}
				}
			}
		}
	}

	// factor CA mask out of CB mask; invert both
	for ( int i=0 ; i<oversamplegrid_[0]*oversamplegrid_[1]*oversamplegrid_[2]; ++i ) {
		rho_cb[i] = (1-rho_cb[i])*rho_ca[i];
		rho_ca[i] = (1-rho_ca[i]);
	}
}


// resample maps subject to rotation
// get self clashes and self rotations
core::Real
CrystDock::resample_maps_and_get_self(
	FArray3D<Real> const &rho_ca, FArray3D<Real> const &rho_cb,
	numeric::xyzMatrix<Real> R, Spacegroup const &/*sg*/, //Commented out by V. Mulligan to fix the debug-mode build.  Sorry to twiddle with your pilot apps, Frank.
	FArray3D<Real> &r_rho_ca, FArray3D<Real> &r_rho_cb,
	utility::vector1<SingleInterface> &p1_interface_map ) {
	r_rho_ca.dimension( grid_[0], grid_[1], grid_[2] ); r_rho_ca=0;
	r_rho_cb.dimension( grid_[0], grid_[1], grid_[2] ); r_rho_cb=0;

	FArray3D<Real> r_rho_ca_base = r_rho_ca;
	FArray3D<Real> r_rho_cb_base = r_rho_cb;

	numeric::xyzMatrix<Real> Ri = numeric::inverse(R);
	numeric::xyzMatrix<Real> Rigridspace = c2i_*Ri*i2c_;
	numeric::xyzMatrix<Real> Rgridspace = c2i_*R*i2c_;

	// rotate "oversmaple grid and find the boundaries
	Real xmax=0, ymax=0, zmax=0;
	{
		numeric::xyzVector<Real> boundbox;
		boundbox = Rgridspace * numeric::xyzVector<Real>(0, 0, oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(0, oversamplegrid_[1], 0);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], 0, 0);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(0, oversamplegrid_[1], oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], 0, oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], oversamplegrid_[1], 0);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
		boundbox = Rgridspace * numeric::xyzVector<Real>(oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2]);
		xmax = std::max( xmax, std::abs(boundbox[0] )); ymax = std::max( ymax, std::abs(boundbox[1] )); zmax = std::max( zmax, std::abs(boundbox[2] ));
	}

	int AMAX = (int)std::ceil( 0.5*(xmax/grid_[0]-1) );
	int BMAX = (int)std::ceil( 0.5*(ymax/grid_[1]-1) );
	int CMAX = (int)std::ceil( 0.5*(zmax/grid_[2]-1) );

	// calculate base transformation
	for ( int z=1; z<=grid_[2]; ++z ) {
		for ( int y=1; y<=grid_[1]; ++y ) {
			for ( int x=1; x<=grid_[0]; ++x ) {
				int cx=x-1,cy=y-1,cz=z-1;
				if ( cx>grid_[0]/2 ) cx -= grid_[0];
				if ( cy>grid_[1]/2 ) cy -= grid_[1];
				if ( cz>grid_[2]/2 ) cz -= grid_[2];

				numeric::xyzVector<Real> rx = Rigridspace*numeric::xyzVector<Real>(cx,cy,cz);

				r_rho_ca_base(x,y,z) = interp_linear( rho_ca, rx );
				r_rho_cb_base(x,y,z) = interp_linear( rho_cb, rx );
			}
		}
	}

	r_rho_ca = r_rho_ca_base;
	r_rho_cb = r_rho_cb_base;

	Real ca_overlap = 0;
	//Real Npoints = grid_[0] * grid_[1] * grid_[2]; //Commented out by V. Mulligan to fix the debug-mode build.  Sorry to twiddle with your pilot apps, Frank.

	// offset transformations
	for ( int a=-AMAX; a<=AMAX; ++a ) {
		for ( int b=-BMAX; b<=BMAX; ++b ) {
			for ( int c=-CMAX; c<=CMAX; ++c ) {
				if ( a==0 && b==0 && c==0 ) continue;

				Real cb_overlap_abc = 0;

				for ( int z=1; z<=grid_[2]; ++z ) {
					for ( int y=1; y<=grid_[1]; ++y ) {
						for ( int x=1; x<=grid_[0]; ++x ) {
							int cx=x-1,cy=y-1,cz=z-1;
							if ( cx>grid_[0]/2 ) cx -= grid_[0];
							if ( cy>grid_[1]/2 ) cy -= grid_[1];
							if ( cz>grid_[2]/2 ) cz -= grid_[2];

							cx += a*grid_[0]; cy += b*grid_[1]; cz += c*grid_[2];

							numeric::xyzVector<Real> rx = Rigridspace*numeric::xyzVector<Real>(cx,cy,cz);

							Real rho_ca_rx = interp_linear( rho_ca, rx );
							Real rho_cb_rx = interp_linear( rho_cb, rx );

							// compute exact overlap
							ca_overlap     += rho_ca_rx*r_rho_ca_base(x,y,z);
							cb_overlap_abc += rho_cb_rx*r_rho_cb_base(x,y,z);

							// add to unit cell (will be used if this is non-overlapping)
							r_rho_ca(x,y,z) += rho_ca_rx;
							r_rho_cb(x,y,z) += rho_cb_rx;
						}
					}
				}

				// add interface if large enough
				cb_overlap_abc *= voxel_volume_;

				if ( cb_overlap_abc > mininterface_ ) {
					p1_interface_map.push_back(
						SingleInterface(
						numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1),
						numeric::xyzVector<core::Real>(a,b,c),
						cb_overlap_abc
						)
					);
				}
			}
		}
	}

	ca_overlap *= voxel_volume_;
	return ca_overlap;
}


Real
CrystDock::get_radius_of_pose( Pose & pose ) {
	Real radius = pose.residue(1).atom(2).xyz().length();
	for ( int i=2; i<=(int)pose.total_residue(); ++i ) {
		Real r=pose.residue(i).atom(2).xyz().length();
		if ( radius < r ) radius=r;
	}
	return radius;
}

Real
CrystDock::get_transform_distance (InterfaceHit ih_vec, InterfaceHit ih_vec_clustered , UniformRotationSampler const &urs, Real radius){
	Real x1=ih_vec.x, y1=ih_vec.y, z1=ih_vec.z, x2=ih_vec_clustered.x, y2=ih_vec_clustered.y, z2=ih_vec_clustered.z;
	Size rot1=ih_vec.rot_index, rot2=ih_vec_clustered.rot_index;

	numeric::xyzMatrix<Real> R1, R2;
	urs.get(rot1,R1);
	urs.get(rot2,R2);

	numeric::xyzVector< core::Real > vec1;
	numeric::xyzVector< core::Real > vec2;

	vec1[0]=x1;vec1[1]=y1;vec1[2]=z1;
	vec2[0]=x2;vec2[1]=y2;vec2[2]=z2;
	Real dist=(vec1-vec2).length();
	R2=numeric::inverse(R2);
	Real angle=DEG2RAD*R2ang(R1*R2);
	Real transform_dist=angle*radius+dist;

	return transform_dist;
}


// TO DO: also align on principal axes to make bounding box as small as possible
numeric::xyzVector<Real>
CrystDock::center_pose_at_origin( Pose & pose ) {
	numeric::xyzVector<Real> com(0,0,0);
	int count=0;
	for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		com += pose.residue(i).atom(2).xyz();
		count++;
	}
	com /= count;

	for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
		for ( int j=1; j<=(int)pose.residue(i).natoms(); ++j ) {
			pose.set_xyz(id::AtomID(j,i), pose.residue(i).atom(j).xyz()-com );
		}
	}
	return com;
}

void
CrystDock::add_crystinfo_to_pose( Pose & pose ) {
	CrystInfo ci;

	ci.A(sg_.A()); ci.B(sg_.B()); ci.C(sg_.C()); ci.alpha(sg_.alpha()); ci.beta(sg_.beta()); ci.gamma(sg_.gamma());
	ci.spacegroup(sg_.pdbname());

	pose.pdb_info()->set_crystinfo(ci);
}

void
CrystDock::dump_transformed_pdb( Pose pose, InterfaceHit ih, UniformRotationSampler const &urs, std::string outname ,std::string basename ) {
	numeric::xyzMatrix<Real> R;
	urs.get( ih.rot_index, R );
	numeric::xyzVector<Real> T (ih.x, ih.y, ih.z);

	// add score to header
	core::pose::RemarkInfo remark;
	std::ostringstream oss;
	oss << "  score = " << ih.score;
	remark.num = 1; remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );

	utility::vector1<std::string> perintinfo = ih.to_string();
	for ( int i=1; i<=(int)perintinfo.size(); ++i ) {
		remark.value = perintinfo[i];
		pose.pdb_info()->remarks().push_back( remark );
	}

	if ( compact_ ) {
		std::ofstream matrix;
		std::string matrixname = basename + std::string(".matrix");
		matrix.open(matrixname.c_str(), std::ios::out | std::ios::app);
		matrix << outname<<"\n";
		matrix << "  score = " << ih.score << "\n";
		matrix << "Rotation matrix: \n";
		matrix << R.xx() << " " << R.xy() << " " << R.xz() << "\n";
		matrix << R.yx() << " " << R.yy() << " " << R.yz() << "\n";
		matrix << R.zx() << " " << R.zy() << " " << R.zz() << "\n";
		matrix << "Transformation Matrix: \n";
		matrix << T[0] << " "<< T[1]<<" "<<T[2]<<std::endl;
		matrix.close();
	} else {
		pose.apply_transform_Rx_plus_v( R,T );
		pose.dump_pdb( outname );
	}

}

// nearest-neighbor interpolation subject to grid-space transform
void
CrystDock::transform_map(
	FArray3D<Real> const &rho,
	numeric::xyzMatrix<Real> S, numeric::xyzVector<Real> T,
	FArray3D<Real> &Srho) {
	Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
	for ( int z=1; z<=rho.u3(); ++z ) {
		for ( int y=1; y<=rho.u2(); ++y ) {
			for ( int x=1; x<=rho.u1(); ++x ) {
				int cx=x-1,cy=y-1,cz=z-1;
				int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) + (T[0]*rho.u1()) + 0.5 ) , rho.u1());
				int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) + (T[1]*rho.u2()) + 0.5 ) , rho.u2());
				int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) + (T[2]*rho.u3()) + 0.5 ) , rho.u3());
				Real rho_sx = rho(rx,ry,rz);
				Srho(x,y,z) = rho_sx;
			}
		}
	}
}

// same as previous function, applies an offset of 0
void
CrystDock::transform_map_offset0(
	FArray3D<Real> const &rho,
	numeric::xyzMatrix<Real> S,
	FArray3D<Real> &Srho) {
	Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
	for ( int z=1; z<=rho.u3(); ++z ) {
		for ( int y=1; y<=rho.u2(); ++y ) {
			for ( int x=1; x<=rho.u1(); ++x ) {
				int cx=x-1,cy=y-1,cz=z-1;
				int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) ) , rho.u1());
				int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) ) , rho.u2());
				int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) ) , rho.u3());
				Real rho_sx = rho(rx,ry,rz);
				Srho(x,y,z) = rho_sx;
			}
		}
	}
}

// wrapper does a series of 2d ffts
//    this could be an option to the low-level fft code to possibly save time
void
CrystDock::fft2dslice( FArray3D<Real> const &rho, FArray3D< std::complex<Real> > &Frho, int axisSlice ) {
	FArray2D<Real> rhoSlice;
	FArray2D< std::complex<Real> > FrhoSlice;
	int xi = rho.u1(), yi = rho.u2(), zi = rho.u3();
	Frho.dimension(xi,yi,zi);

	if ( axisSlice == 1 ) {
		rhoSlice.dimension(yi,zi);
		for ( int ii=1; ii<=xi; ii++ ) {
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rhoSlice(jj,kk) = rho(ii,jj,kk);
				}
			}
			numeric::fourier::fft2(rhoSlice, FrhoSlice);
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					Frho(ii,jj,kk) = FrhoSlice(jj,kk);
				}
			}
		}
	} else if ( axisSlice == 2 ) {
		rhoSlice.dimension(xi,zi);
		for ( int jj=1; jj<=yi; jj++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rhoSlice(ii,kk) = rho(ii,jj,kk);
				}
			}
			numeric::fourier::fft2(rhoSlice, FrhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					Frho(ii,jj,kk) = FrhoSlice(ii,kk);
				}
			}
		}
	} else if ( axisSlice == 3 ) {
		rhoSlice.dimension(xi,yi);
		for ( int kk=1; kk<=zi; kk++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					rhoSlice(ii,jj) = rho(ii,jj,kk);
				}
			}
			numeric::fourier::fft2(rhoSlice, FrhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					Frho(ii,jj,kk) = FrhoSlice(ii,jj);
				}
			}
		}
	} else {
		utility_exit_with_message( "ERROR! Bad axis specified!");
	}
}

// wrapper does a series of 2d iffts
void
CrystDock::ifft2dslice( FArray3D< std::complex<Real> > const &Frho, FArray3D<Real> &rho, int axisSlice ) {
	FArray2D<Real> rhoSlice;
	FArray2D< std::complex<Real> > FrhoSlice;
	int xi = Frho.u1(), yi = Frho.u2(), zi = Frho.u3();
	rho.dimension(xi,yi,zi);
	rho = 0;

	if ( axisSlice == 1 ) {
		FrhoSlice.dimension(yi,zi);
		for ( int ii=1; ii<=xi; ii++ ) {
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					FrhoSlice(jj,kk) = Frho(ii,jj,kk);
				}
			}
			numeric::fourier::ifft2(FrhoSlice, rhoSlice);
			for ( int jj=1; jj<=yi; jj++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rho(ii,jj,kk) = rhoSlice(jj,kk);
				}
			}
		}
	} else if ( axisSlice == 2 ) {
		FrhoSlice.dimension(xi,zi);
		for ( int jj=1; jj<=yi; jj++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					FrhoSlice(ii,kk) = Frho(ii,jj,kk);
				}
			}
			numeric::fourier::ifft2(FrhoSlice, rhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int kk=1; kk<=zi; kk++ ) {
					rho(ii,jj,kk) = rhoSlice(ii,kk);
				}
			}
		}
	} else if ( axisSlice == 3 ) {
		FrhoSlice.dimension(xi,yi);
		for ( int kk=1; kk<=zi; kk++ ) {
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					FrhoSlice(ii,jj) = Frho(ii,jj,kk);
				}
			}
			numeric::fourier::ifft2(FrhoSlice, rhoSlice);
			for ( int ii=1; ii<=xi; ii++ ) {
				for ( int jj=1; jj<=yi; jj++ ) {
					rho(ii,jj,kk) = rhoSlice(ii,jj);
				}
			}
		}
	} else {
		utility_exit_with_message( "ERROR! Bad axis specified!");
	}
}

void
CrystDock::project_along_axis( FArray3D<Real> &rho, numeric::xyzVector<Real> axis ) {
	FArray3D<Real> rho_input = rho;

	Real scalefact = std::max(std::fabs(axis[0]), std::max(std::fabs(axis[1]),std::fabs(axis[2])));
	axis /= scalefact;

	int ngrid = (std::fabs(axis[0]>0.999))? grid_[0] : ((std::fabs(axis[1])>0.999)? grid_[1] : grid_[2]);

	// for pp = 1:D-1; mfft = mfft+circshift( mfft_orig, [pp*slide(1,1),pp*sl ide(2,2),pp*slide(3,3)]);
	for ( int P=1; P<ngrid; ++P ) {  // one less than full cycle
		for ( int z=1; z<=grid_[2]; ++z ) {
			for ( int y=1; y<=grid_[1]; ++y ) {
				for ( int x=1; x<=grid_[0]; ++x ) {
					int xt = pos_mod((int)floor( x+P*axis[0]-0.5 ), grid_[0]) + 1;
					int yt = pos_mod((int)floor( y+P*axis[1]-0.5 ), grid_[1]) + 1;
					int zt = pos_mod((int)floor( z+P*axis[2]-0.5 ), grid_[2]) + 1;
					rho(x,y,z) += rho_input(xt,yt,zt);
				}
			}
		}
	}
}

// the main convolution operation used for both CA and CB maps
void
CrystDock::do_convolution( FArray3D<Real> const &rho, FArray3D<Real> const &Srho, numeric::xyzMatrix<Real> S, FArray3D<Real> &conv_out) {
	FArray3D<Real> working_s, working_trans, working_strans;
	FArray3D< std::complex<Real> > Fworking_trans, Fworking_strans;
	Size Npoints = rho.u1()*rho.u2()*rho.u3();

	// find out what plane we are slicing along
	Quat Sinv_Q(S);

	// many different possibilities explain the convoluted logic below
	//    - this may return (1,0,0), (sqrt(1/2),sqrt(1/2),0) or (1/2,1/2,1/2) for X, XY, and XYZ symm axes
	//    - also in some spacegps( e.g. P3112[5] ) it may be (1,1/4,0) -- since it's not a rotation but a skew
	//    - in others ( e.g. P6122[8] ) it may also be (1,1/4,0) but we need to pre-skew
	//    - in these cases, we want to process as 'x' rotation with the skew handled in the projection
	int slice_x = (Sinv_Q.x_>=0.5) ? 1 : ((Sinv_Q.x_<=-0.5) ? -1:0);
	int slice_y = (Sinv_Q.y_>=0.5) ? 1 : ((Sinv_Q.y_<=-0.5) ? -1:0);
	int slice_z = (Sinv_Q.z_>=0.5) ? 1 : ((Sinv_Q.z_<=-0.5) ? -1:0);

	numeric::xyzVector<Real> rotaxis(
		std::sqrt(std::fabs(Sinv_Q.x_)),
		std::sqrt(std::fabs(Sinv_Q.y_)),
		std::sqrt(std::fabs(Sinv_Q.z_)));
	if ( Sinv_Q.x_<0 ) rotaxis[0] *= -1;
	if ( Sinv_Q.y_<0 ) rotaxis[1] *= -1;
	if ( Sinv_Q.z_<0 ) rotaxis[2] *= -1;

	// do we need to do fft convolution at all?
	if ( slice_x==0 && slice_y==0 && slice_z==0 ) {

		// translation only, just take dot product

		if ( debug_||debug_exact_ ) {
			TR << "no slice"  << std::endl;
			TR << "                : " << Sinv_Q.x_ << " " << Sinv_Q.y_ << " " << Sinv_Q.z_ << std::endl;
		}
		core::Real dotProd=0;
		for ( int i=0; i<(int)Npoints; ++i ) dotProd += rho[i]*Srho[i];
		conv_out.dimension( rho.u1(),rho.u2(),rho.u3() );
		conv_out = dotProd;
	} else {

		// do the fft convolution

		// Rinv maps xyz->pij (our reindexed coordinates with p is perpendicular to the symmaxis)
		// R maps pij->xyz
		numeric::xyzMatrix<Real> Rinv;
		int slice_axis = 0;
		bool R_is_identity=false;
		if ( slice_x != 0 ) {
			slice_axis = 1;
			Rinv.row_x( Vector( slice_x, slice_y, slice_z ) );
			R_is_identity = (slice_y==0 && slice_z==0);
		} else {
			Rinv.row_x( Vector( 1, 0, 0 ) );
		}

		if ( slice_x == 0 && slice_y != 0 ) {
			slice_axis = 2;
			Rinv.row_y( Vector( slice_x, slice_y, slice_z ) );
			R_is_identity = (slice_z==0);
		} else {
			Rinv.row_y( Vector( 0, 1, 0 ) );
		}
		if ( slice_x == 0 && slice_y == 0 && slice_z != 0 ) {
			slice_axis = 3;
			Rinv.row_z( Vector( slice_x, slice_y, slice_z ) );
			R_is_identity = true;
		} else {
			Rinv.row_z( Vector( 0, 0, 1 ) );
		}

		// oddball case (+ rotations of it)
		if ( S == numeric::xyzMatrix<Real>::rows(1,-1,0,   0,-1,0,   0,0,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 2,-1,0,  1,0,0,  0,0,1 );
			R_is_identity = false;
			rotaxis = Vector(1,0,0);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(-1,0,0, -1,1,0, 0,0,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, -1,2,0, 0,0,1 );
			R_is_identity = false;
			rotaxis = Vector(0,1,0);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(-1,0,0, 0,-1,0, 0,-1,1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows(1,0,0, 0,0,1, 0,-1,2);
			R_is_identity = false;
			rotaxis = Vector(0,0,1);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(1,0,-1, 0,-1,0, 0,0,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 2,0,-1, 0,1,0, 1,0,0 );
			R_is_identity = false;
			rotaxis = Vector(1,0,0);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(0,0,-1, -1,0,0, 0,1,-1) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, 1,0,0, 0,2,-1 );
			R_is_identity = false;
			rotaxis = Vector(0,0,1);
		}
		if ( S == numeric::xyzMatrix<Real>::rows(0,-1,0, -1,0,1, -1,0,0) ) {
			Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, -1,0,2, 0,0,1 );
			R_is_identity = false;
			rotaxis = Vector(0,1,0);
		}

		numeric::xyzMatrix<Real> R = numeric::inverse(Rinv);

		// Q is our symmetric transformation in our new coordinates
		numeric::xyzMatrix<Real> Q = Rinv*S*R;

		// Qstar is the FFT transformation
		//   Q-I on diagonal elements not on the slice axis
		//   partially apply the elements in the column of the slice axis (only applied to symm-sampled map)
		//        to compensate for point group sliding as we move up the symm axis in resampled groups
		numeric::xyzMatrix<Real> QstarF, QstarG;
		QstarF.clear(); //SML 04-25-14; my compiler warns Qstar* may be used uninitialized (because all inits below are in ifs)
		QstarG.clear();
		if ( slice_axis == 1 ) {
			QstarF = numeric::xyzMatrix<Real>::rows( 1,0,0,   Q.yx(),Q.yy()-1,Q.yz(),   Q.zx(),Q.zy(),Q.zz()-1 );
			QstarG = numeric::xyzMatrix<Real>::rows( 1,0,0,        0,Q.yy()-1,Q.yz(),        0,Q.zy(),Q.zz()-1 );
		} else if ( slice_axis == 2 ) {
			QstarF = numeric::xyzMatrix<Real>::rows( Q.xx()-1,  Q.xy(),Q.xz(), 0,1,0, Q.zx(),  Q.zy(),Q.zz()-1 );
			QstarG = numeric::xyzMatrix<Real>::rows( Q.xx()-1,       0,Q.xz(), 0,1,0, Q.zx(),       0,Q.zz()-1 );
		} else if ( slice_axis == 3 ) {
			QstarF = numeric::xyzMatrix<Real>::rows( Q.xx()-1,Q.xy(),  Q.xz(), Q.yx(),Q.yy()-1,  Q.yz(), 0,0,1 );
			QstarG = numeric::xyzMatrix<Real>::rows( Q.xx()-1,Q.xy(),       0, Q.yx(),Q.yy()-1,       0, 0,0,1 );
		}

		if ( debug_||debug_exact_ ) {
			TR << "  slice on axis :    " << slice_axis << std::endl;
			TR << "slice along axis: " << slice_x << " " << slice_y << " " << slice_z << std::endl;
			TR << "                : " << Sinv_Q.x_ << " " << Sinv_Q.y_ << " " << Sinv_Q.z_ << std::endl;
			TR << "                : " << rotaxis[0] << " " << rotaxis[1] << " " << rotaxis[2] << std::endl;
			TR << "         S = [" << S.xx() << "," <<  S.xy() << "," <<  S.xz() << ";  "
				<< S.yx() << "," <<  S.yy() << "," <<  S.yz() << ";  "
				<< S.zx() << "," <<  S.zy() << "," <<  S.zz() << "]" << std::endl;
			if ( !R_is_identity ) {
				TR << "      Rinv = [" << Rinv.xx() << "," << Rinv.xy() << "," << Rinv.xz() << ";  "
					<< Rinv.yx() << "," << Rinv.yy() << "," << Rinv.yz() << ";  "
					<< Rinv.zx() << "," << Rinv.zy() << "," << Rinv.zz() << "]" << std::endl;
				TR << "         R = [" << R.xx() << "," <<  R.xy() << "," <<  R.xz() << ";  "
					<< R.yx() << "," <<  R.yy() << "," <<  R.yz() << ";  "
					<< R.zx() << "," <<  R.zy() << "," <<  R.zz() << "]" << std::endl;
			}
			TR << "         Q = [" << Q.xx() << "," << Q.xy() << "," << Q.xz() << ";  "
				<< Q.yx() << "," << Q.yy() << "," << Q.yz() << ";  "
				<< Q.zx() << "," << Q.zy() << "," << Q.zz() << "]" << std::endl;
			TR << "        Qf = [" << QstarF.xx() << "," <<  QstarF.xy() << "," <<  QstarF.xz() << ";  "
				<< QstarF.yx() << "," <<  QstarF.yy() << "," <<  QstarF.yz() << ";  "
				<< QstarF.zx() << "," <<  QstarF.zy() << "," <<  QstarF.zz() << "]" << std::endl;
			TR << "        Qg = [" << QstarG.xx() << "," <<  QstarG.xy() << "," <<  QstarG.xz() << ";  "
				<< QstarG.yx() << "," <<  QstarG.yy() << "," <<  QstarG.yz() << ";  "
				<< QstarG.zx() << "," <<  QstarG.zy() << "," <<  QstarG.zz() << "]" << std::endl;
		}

		// transform
		transform_map_offset0( rho, R, working_trans);
		transform_map_offset0( Srho, R, working_strans);

		// transform in each plane
		// fft in each plane
		FArray3D<Real> rworking_trans, rworking_strans;
		transform_map_offset0( working_trans, QstarF, rworking_trans);
		transform_map_offset0( working_strans, QstarG, rworking_strans);

		fft2dslice( rworking_trans, Fworking_trans, slice_axis );
		fft2dslice( rworking_strans, Fworking_strans, slice_axis );
		for ( int i=0; i<(int)Npoints; ++i ) Fworking_trans[i] *= std::conj(Fworking_strans[i]);


		if ( R_is_identity ) {
			ifft2dslice(Fworking_trans, conv_out, slice_axis);
		} else {
			ifft2dslice(Fworking_trans, working_trans, slice_axis);
			transform_map_offset0( working_trans, Rinv, conv_out);
		}
		project_along_axis( conv_out, rotaxis );
	}
}


// given a list of symmops that are contacting the main subunit,
//    break down the contribution into individual interfaces
//    find the weakest connection in fully connected lattice
// work in the space of oversample_grid_
Real
CrystDock::get_interfaces(
	utility::vector1<core::kinematics::RT> const &rts,
	numeric::xyzMatrix<Real> R,
	numeric::xyzVector<Real> xyz_grid,
	utility::vector1<Size> symmopList,
	FArray3D<Real> const &rho_cb,
	utility::vector1<SingleInterface> &allInterfaces ) {

	numeric::xyzMatrix<Real> Rgridspace = c2i_*R*i2c_;
	numeric::xyzMatrix<Real> R_inv = numeric::inverse(Rgridspace);

	// [STAGE 1] compute exact interface scores
	Real best_int = 0;
	Real cb_sum = 0;
	for ( int i=1; i<=(int)symmopList.size(); ++i ) {
		int s = symmopList[i];

		numeric::xyzMatrix<Real> s_i = rts[s].get_rotation();
		numeric::xyzVector<Real> t_i = rts[s].get_translation();
		numeric::xyzMatrix<Real> s_inv = numeric::inverse(s_i);

		// find the symm offset to minimize x->Sx distance
		int A0=0,B0=0,C0=0;
		Real mindist = 1e30;
		for ( int a=-3; a<=3; ++a ) {
			for ( int b=-3; b<=3; ++b ) {
				for ( int c=-3; c<=3; ++c ) {
					numeric::xyzVector<Real> A((t_i[0]+a)*grid_[0],(t_i[1]+b)*grid_[1],(t_i[2]+c)*grid_[2]);
					Real dist_abc = (xyz_grid - (s_i*(xyz_grid) + A)).length_squared();
					if ( dist_abc < mindist ) {
						mindist = dist_abc;
						A0=a;B0=b;C0=c;
					}
				}
			}
		}

		//Real Npoints = grid_[0]*grid_[1]*grid_[2]; //Commented out by V. Mulligan to fix the debug-mode build.  Sorry to twiddle with your pilot apps, Frank.

		// now look +/-1 copy in each direction
		for ( int a=A0-1; a<=A0+1; ++a ) {
			for ( int b=B0-1; b<=B0+1; ++b ) {
				for ( int c=C0-1; c<=C0+1; ++c ) {
					Real cb_overlap_abc = 0;
					for ( int z=1; z<=oversamplegrid_[2]; ++z ) {
						for ( int y=1; y<=oversamplegrid_[1]; ++y ) {
							for ( int x=1; x<=oversamplegrid_[0]; ++x ) {
								int cx=x-1,cy=y-1,cz=z-1;
								if ( cx>oversamplegrid_[0]/2 ) cx -= oversamplegrid_[0];
								if ( cy>oversamplegrid_[1]/2 ) cy -= oversamplegrid_[1];
								if ( cz>oversamplegrid_[2]/2 ) cz -= oversamplegrid_[2];

								// given:
								//    R,C = original rotation & translation
								//    S,T,A = symmetric transform & unit cell offset
								// we want overlap of (Rx+C) and (S*(Rx+C)+T+A)
								// ==> offset of x and R^-1*( S*(Rx+C)+T+A - C )
								numeric::xyzVector<Real> x_i(cx,cy,cz);
								numeric::xyzVector<Real> A((t_i[0]+a)*grid_[0],(t_i[1]+b)*grid_[1],(t_i[2]+c)*grid_[2]);
								numeric::xyzVector<Real> transformX = R_inv*(s_i*(Rgridspace*x_i + xyz_grid) + A - xyz_grid);

								cb_overlap_abc += rho_cb(x,y,z) * interp_linear( rho_cb, transformX );
							}
						}
					}

					// add interface if large enough
					cb_overlap_abc *= voxel_volume_;
					cb_sum += cb_overlap_abc;
					best_int = std::max( best_int, cb_overlap_abc );
					if ( cb_overlap_abc > mininterface_ ) {
						allInterfaces.push_back(
							SingleInterface( s_i, numeric::xyzVector<core::Real>(t_i[0]+a,t_i[1]+b,t_i[2]+c), cb_overlap_abc ) );
					}
				}
			}
		}
	}
	return cb_sum;
}

// get the score per interface
void
CrystDock::get_interfaces_allatom(
	Pose pose, // make a copy
	utility::vector1<core::kinematics::RT> const &/*rts*/,
	numeric::xyzMatrix<Real> R,
	numeric::xyzVector<Real> xyz,
	utility::vector1<SingleInterface> &allInterfaces )
{
	pose.apply_transform_Rx_plus_v( R, xyz );

	//////////////////////
	////
	//// recenter
	////
	//////////////////////
	Size nres = pose.total_residue();

	Vector com(0,0,0);
	for ( Size i=1; i<= nres; ++i ) com += pose.residue(i).xyz(2);
	com /= nres;

	Real mindis2(1e6);
	for ( Size i=1; i<= nres; ++i ) {
		Real const dis2( com.distance_squared(  pose.residue(i).xyz(2) ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
		}
	}
	Size nsymm = sg_.nsymmops();
	Size bestxform=0;
	Vector bestoffset(0,0,0);
	mindis2=1e6;
	com = sg_.c2f()*com;
	for ( Size i=1; i<=nsymm; ++i ) {
		Vector foffset = sg_.symmop(i).get_rotation()*com + sg_.symmop(i).get_translation(), rfoffset;
		rfoffset[0] = min_mod( foffset[0], 1.0 );
		rfoffset[1] = min_mod( foffset[1], 1.0 );
		rfoffset[2] = min_mod( foffset[2], 1.0 );
		Real dist = (sg_.f2c()*rfoffset).length_squared();
		if ( dist<mindis2 ) {
			mindis2=dist;
			bestxform=i;
			bestoffset = foffset - rfoffset;
		}
	}
	numeric::xyzMatrix<Real> Rx = sg_.f2c()*sg_.symmop(bestxform).get_rotation()*sg_.c2f();
	numeric::xyzVector<Real> Tx = sg_.f2c()*(sg_.symmop(bestxform).get_translation() - bestoffset);
	pose.apply_transform_Rx_plus_v( Rx,Tx );

	//////////////////////
	////
	//// find contacts
	////
	//////////////////////
	Real contact_dist=12;
	Real radius = 0;
	utility::vector1<Vector> monomer_cbs;

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		Size atm=(pose.residue(i).aa() == core::chemical::aa_gly)?2:5;
		Vector cb_i = pose.residue(i).xyz(atm);
		monomer_cbs.push_back(cb_i);
		radius = std::max( (cb_i).length_squared() , radius );
	}
	Size nres_monomer = monomer_cbs.size();
	radius = sqrt(radius);

	for ( int s=1; s<=(int)sg_.nsymmops(); ++s ) {
		numeric::xyzMatrix<Real> R_i = sg_.symmop(s).get_rotation();

		for ( int i=-1; i<=1; ++i ) {
			for ( int j=-1; j<=1; ++j ) {
				for ( int k=-1; k<=1; ++k ) {
					if ( s==1 && i==0 && j==0 && k==0 ) continue;

					numeric::xyzVector<Real> T_i = sg_.symmop(s).get_translation() + numeric::xyzVector<Real>(i,j,k);

					// pass 1 check vrt-vrt dist to throw out very distant things
					Real disVRT = T_i.length();
					if ( disVRT>contact_dist+2*radius ) continue;

					// pass 2 check ca-ca dists
					numeric::xyzMatrix<core::Real> R_i_realspace =  sg_.f2c()*R_i*sg_.c2f();
					Size contact=0;
					for ( Size jj=1; jj<= nres_monomer; ++jj ) {
						Vector Xi = R_i_realspace*monomer_cbs[jj] + sg_.f2c()*T_i;
						for ( Size kk=1; kk<= nres_monomer; ++kk ) {
							if ( (Xi-monomer_cbs[kk]).length_squared() < contact_dist*contact_dist ) contact++;
						}
					}

					if ( contact>mininterface_ ) {
						// make minipose & score
						core::pose::Pose poseCopy = pose;
						allInterfaces.push_back(
							SingleInterface( R_i, T_i, (Real)contact ) );
					}
				}
			}
		}
	}
}


core::Real
CrystDock::get_interface_score(
	Pose pose,
	utility::vector1<SingleInterface> &allInterfaces,
	utility::vector1<core::kinematics::RT> const &rts_all ) {
	// [STAGE 2] see if we form a connected lattice
	//    * idea: expand all contacting symmops
	//    *       ensure we generate all symmops
	//    *       ensure we generate (+/-1,0,0), (0,+/-1,0), and (0,0,+/-1)
	int EXPAND_ROUNDS=option[crystdock::expand_rounds];  // no idea if this is sufficient

	int nxformsOrig = (int)allInterfaces.size();

	// stopping conditions
	if ( debug_ || debug_exact_ ) {
		TR << "[EXPAND] Round 0: have " << allInterfaces.size() << " subunits" << std::endl;
		for ( int i=1; i<=(int)allInterfaces.size(); ++i ) {
			TR << "model " << i<< std::endl;
			TR << "[SCORE=" <<  allInterfaces[i].cb_overlap_ << "] R = ["
				<< allInterfaces[i].R_.xx() << "," << allInterfaces[i].R_.xy() << "," << allInterfaces[i].R_.xz() << ";"
				<< allInterfaces[i].R_.yx() << "," << allInterfaces[i].R_.yy() << "," << allInterfaces[i].R_.yz() << ";"
				<< allInterfaces[i].R_.zx() << "," << allInterfaces[i].R_.zy() << "," << allInterfaces[i].R_.zz()
				<< "]";
			TR << " T = [" << allInterfaces[i].T_[0] << "," << allInterfaces[i].T_[1] << "," << allInterfaces[i].T_[2] << "]" << std::endl;
		}
	}

	for ( int rd=1; rd<=EXPAND_ROUNDS; ++rd ) {
		int nxforms = (int)allInterfaces.size();
		for ( int i=1; i<=nxforms; ++i ) {
			for ( int j=1; j<=nxformsOrig; ++j ) {
				// keep in inner loop
				numeric::xyzMatrix<Real> const &Si=allInterfaces[i].R_;
				numeric::xyzVector<Real> const &Ti=allInterfaces[i].T_;
				numeric::xyzMatrix<Real> const &Sj=allInterfaces[j].R_;
				numeric::xyzVector<Real> const &Tj=allInterfaces[j].T_;

				numeric::xyzMatrix<Real> Sij = Si*Sj;
				numeric::xyzVector<Real> Tij = Ti + Si*Tj;
				Real overlap_ij = std::min( allInterfaces[i].cb_overlap_, allInterfaces[j].cb_overlap_ );

				if ( (Tij[0])>2 || (Tij[1])>2 || (Tij[2])>2 ) continue;
				if ( (Tij[0])<-2 || (Tij[1])<-2 || (Tij[2])<-2 ) continue;

				// check if it is in the set
				// we could hash these for a small speed increase...
				bool unique_xform = true;
				for ( int k=1; k<=(int)allInterfaces.size(); ++k ) {
					if ( transforms_equiv(Sij, Tij, allInterfaces[k].R_, allInterfaces[k].T_) ) {
						// check if min interface is better
						allInterfaces[k].cb_overlap_ = std::max( allInterfaces[k].cb_overlap_, overlap_ij );
						unique_xform = false;
						//TR << "k= " << k << ", allInterfaces.size= " << (int)allInterfaces.size() << std::endl;
					}
				}
				if ( unique_xform ) {
					allInterfaces.push_back( SingleInterface( Sij, Tij, overlap_ij ) );
				}
			}
		}

		if ( debug_ || debug_exact_ ) {
			if ( rd==1 ) {
				TR << "[EXPAND] Round " << rd << ": have " << allInterfaces.size() << " subunits" << std::endl;
				for ( int i=1; i<=(int)allInterfaces.size(); ++i ) {
					TR << "[SCORE=" <<  allInterfaces[i].cb_overlap_ << "] R = ["
						<< allInterfaces[i].R_.xx() << "," << allInterfaces[i].R_.xy() << "," << allInterfaces[i].R_.xz() << ";"
						<< allInterfaces[i].R_.yx() << "," << allInterfaces[i].R_.yy() << "," << allInterfaces[i].R_.yz() << ";"
						<< allInterfaces[i].R_.zx() << "," << allInterfaces[i].R_.zy() << "," << allInterfaces[i].R_.zz()
						<< "]";
					TR << " T = [" << allInterfaces[i].T_[0] << "," << allInterfaces[i].T_[1] << "," << allInterfaces[i].T_[2] << "]" << std::endl;

					Pose poseCopy = pose;
					numeric::xyzVector<Real> Tnew = i2c_*numeric::xyzVector<Real>(
						(Real)grid_[0]*allInterfaces[i].T_[0],
						(Real)grid_[1]*allInterfaces[i].T_[1],
						(Real)grid_[2]*allInterfaces[i].T_[2]);

					poseCopy.apply_transform_Rx_plus_v( i2c_*allInterfaces[i].R_*c2i_,Tnew);
					poseCopy.append_pose_by_jump (pose, 1);
					std::ostringstream oss;
					oss << "model_" << i << ".pdb";
					poseCopy.dump_pdb( oss.str() );
				}
			}
		}
	}


	// check if we are fully connected
	//   1--- all symmops
	Real score = 1e30;
	bool connected = true;
	for ( int i=2; i<=(int)rts_all.size() && connected; ++i ) {  // 1 is identity
		bool contains_j = false;
		for ( int j=1; j<=(int)allInterfaces.size() && !contains_j; ++j ) {
			numeric::xyzMatrix<Real> const &Si=allInterfaces[j].R_;
			numeric::xyzVector<Real> const &Ti=allInterfaces[j].T_;
			if ( transforms_equiv( Si, Ti, rts_all[i].get_rotation(), rts_all[i].get_translation()) ) {
				contains_j = true;
				score = std::min( score, allInterfaces[j].cb_overlap_ );
			}
		}
		//TR << "Connection check  at " << i << "and " << contains_j <<std::endl;
		connected &= contains_j;
	}
	if ( !connected ) return 0.0;

	//   2--- +/- 1 in each direction
	Real score_00p=0, score_00m=0, score_0m0=0, score_0p0=0, score_p00=0, score_m00=0;
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	for ( int i=1; i<=(int)allInterfaces.size(); ++i ) {
		numeric::xyzMatrix<Real> const &Si=allInterfaces[i].R_;
		numeric::xyzVector<Real> const &Ti=allInterfaces[i].T_;
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0, 1)) ) {
			score_00p = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0,-1)) ) {
			score_00m = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0, 1,0)) ) {
			score_0p0 = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,-1,0)) ) {
			score_0m0 = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>( 1,0,0)) ) {
			score_p00 = allInterfaces[i].cb_overlap_;
		}
		if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(-1,0,0)) ) {
			score_m00 = allInterfaces[i].cb_overlap_;
		}
	}
	//TR << "Scores: score_00p " << score_00p << " score_00m: " << score_00m << " score_0m0: " << score_0m0 << " score_0p0: " << score_0p0 << " score_p00: " << score_p00 << " score_m00: " << score_m00 <<std::endl;
	score = std::min( score, std::min(score_00p, score_00m));
	score = std::min( score, std::min(score_0m0, score_0p0));
	score = std::min( score, std::min(score_p00, score_m00));

	return score;
}

core::Real
CrystDock::get_clash_score_exact(
	numeric::xyzVector<int> xyz_grid,
	numeric::xyzMatrix<Real> R,
	numeric::xyzVector<Real> T,
	FArray3D<Real> const &r_rho_ca
) {
	FArray3D<Real> shift_rho_ca, s_shift_rho_ca;
	shift_rho_ca.dimension( grid_[0] , grid_[1] , grid_[2] );
	Size Npoints = grid_[0]*grid_[1]*grid_[2];

	for ( int z=1; z<=grid_[2]; ++z ) {
		for ( int y=1; y<=grid_[1]; ++y ) {
			for ( int x=1; x<=grid_[0]; ++x ) {
				int cx = pos_mod( x-xyz_grid[0]-1 , grid_[0] ) + 1;
				int cy = pos_mod( y-xyz_grid[1]-1 , grid_[1] ) + 1;
				int cz = pos_mod( z-xyz_grid[2]-1 , grid_[2] ) + 1;
				shift_rho_ca(x,y,z) = r_rho_ca(cx,cy,cz);
			}
		}
	}
	transform_map( shift_rho_ca, R,T, s_shift_rho_ca);

	Real retval = 0;
	for ( int i=0; i<(int)Npoints; ++i ) {
		retval += shift_rho_ca[i] * s_shift_rho_ca[i];
	}

	return retval*voxel_volume_;
}


void
CrystDock::apply( Pose & pose) {
	// set up crystal info
	sg_.set_spacegroup( option[ crystdock::spacegroup ] );
	sg_.set_parameters(
		option[ crystdock::A ],option[ crystdock::B ],option[ crystdock::C ],
		option[ crystdock::alpha ],option[ crystdock::beta ],option[ crystdock::gamma ] );

	Real rotstep = option[ crystdock::rot_step ];
	Real maxclash = option[ crystdock::maxclash ];
	// center pose at origin
	numeric::xyzVector<Real> native_shift = center_pose_at_origin( pose );
	add_crystinfo_to_pose( pose );
	int translation_count=1;
	// random reorient
	if ( option[ crystdock::randomize_orientation ]() ) {
		numeric::xyzMatrix<Real> Rrand = protocols::geometry::random_reorientation_matrix( 360.0, 360.0 );
		pose.apply_transform_Rx_plus_v( Rrand, numeric::xyzVector<Real>(0,0,0) );
	} else if ( option[ crystdock::random_rotate ]() > 0 ) {
		const double psi( 2 * numeric::constants::d::pi * numeric::random::uniform() );
		const double theta( std::acos(numeric::sin_cos_range( 1.0 - 2.0  *numeric::random::uniform() ) ) );
		Vector axis( sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta) );
		numeric::xyzMatrix<Real> Rrand = numeric::rotation_matrix( axis, DEG2RAD*option[ crystdock::random_rotate ]() );
		pose.apply_transform_Rx_plus_v( Rrand, numeric::xyzVector<Real>(0,0,0) );
	}


	// get SS
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// lookup symmops
	utility::vector1<core::kinematics::RT> rts = sg_.symmops();
	CheshireCell cc = sg_.cheshire_cell();

	if ( debug_||debug_exact_ ) {
		TR << "search [ " << cc.low[0] <<","<< cc.low[1] <<","<< cc.low[2] << " : "
			<< cc.high[0] <<","<< cc.high[1] <<","<< cc.high[2] << "]" << std::endl;
	}

	// compute rho_ca, rho_cb
	FArray3D<Real> rho_ca, rho_cb, r_rho_ca, r_rho_cb;


	setup_maps( pose, rho_ca, rho_cb, trans_step_);

	Size Npoints = grid_[0]*grid_[1]*grid_[2];
	voxel_volume_ = sg_.volume() / Npoints;

	numeric::xyzVector<Size> ccIndexLow(
		(Size)std::floor(cc.low[0]*grid_[0]+1.5),
		(Size)std::floor(cc.low[1]*grid_[1]+1.5),
		(Size)std::floor(cc.low[2]*grid_[2]+1.5));
	numeric::xyzVector<Size> ccIndexHigh(
		(Size)std::floor(cc.high[0]*grid_[0]+0.5),
		(Size)std::floor(cc.high[1]*grid_[1]+0.5),
		(Size)std::floor(cc.high[2]*grid_[2]+0.5));

	for ( int i=0; i<3; i++ ) {
		if ( ccIndexHigh[i]<ccIndexLow[i] ) {
			ccIndexHigh[i]=ccIndexLow[i];
		}
	}

	if ( debug_||debug_exact_ ) {
		writeMRC( rho_ca, "ca_mask.mrc", true, true );
		writeMRC( rho_cb, "cb_mask.mrc", true, true );
	}

	// space for intermediate results
	FArray3D<Real> working_s, conv_out;

	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	numeric::xyzMatrix<Real> r_local;

	// the collection of hits
	InterfaceHitDatabase IDB(10*nmodels_);
	core::Size nnonclashing = 0, nconnected = 0;

	// foreach rotation
	UniformRotationSampler urs( rotstep );
	urs.remove_redundant( rts, i2c_, c2i_ );

	Size rot_lb = 1, rot_ub = urs.nrots();
	if ( option[ crystdock::rotnum ].user() ) {
		rot_lb = rot_ub = option[ crystdock::rotnum ]();
	}
	Size sym_lb = 2, sym_ub = rts.size();
	if ( option[ crystdock::symnum ].user() ) {
		sym_lb = sym_ub = option[ crystdock::symnum ]();
	}

	// hook into other code for evaluating native interfaces
	if ( eval_native_ ) {
		mininterface_ = 0.0001;  // we want to evaluate all interfaces

		utility::vector1<Size> symmoplist;

		numeric::xyzVector<Real> offset_grid;
		numeric::xyzVector<int> offset_grid_pt;
		utility::vector1<SingleInterface> iinfo;

		offset_grid = c2i_*native_shift;

		// round to nearest grid
		offset_grid_pt[0] = (int) std::floor( offset_grid[0] + 0.5 );
		offset_grid_pt[1] = (int) std::floor( offset_grid[1] + 0.5 );
		offset_grid_pt[2] = (int) std::floor( offset_grid[2] + 0.5 );
		offset_grid = numeric::xyzVector<Real>(offset_grid_pt[0], offset_grid_pt[1], offset_grid_pt[2]);
		native_shift = i2c_*offset_grid;

		// exact ca overlap volume
		core::Real ca_score = 0;
		ca_score += resample_maps_and_get_self( rho_ca, rho_cb, identity, sg_, r_rho_ca, r_rho_cb, iinfo );
		for ( int s=2; s<=(int)rts.size(); ++s ) {
			symmoplist.push_back(s);
			ca_score += get_clash_score_exact( offset_grid_pt, rts[s].get_rotation(), rts[s].get_translation(), r_rho_ca );
		}

		// exact cb overlap sum
		core::Real cb_sum=0;
		for ( int i=1; i<=(int)iinfo.size(); ++i ) {
			cb_sum+= iinfo[i].cb_overlap_;
		}
		cb_sum += get_interfaces( rts, identity, offset_grid, symmoplist, rho_cb, iinfo );

		// exact cb contact count
		iinfo.clear();
		get_interfaces_allatom( pose, rts, identity, native_shift, iinfo );
		core::Real cb_score = get_interface_score( pose, iinfo, rts );

		TR << "OVERLAP score = " << ca_score << std::endl;
		TR << "INTERACTION sum (density) = " << cb_sum << std::endl;
		TR << "INTERACTION score (cb count)= " << cb_score << std::endl;

		// finally shift the pose + write to disk
		pose.apply_transform_Rx_plus_v( identity, i2c_*numeric::xyzVector<Real>(offset_grid_pt[0],offset_grid_pt[1],offset_grid_pt[2]) );

		numeric::xyzVector<Real> xyz = i2c_*numeric::xyzVector<Real>(offset_grid_pt[0],offset_grid_pt[1],offset_grid_pt[2]);

		std::string base_name = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();
		utility::vector1< std::string > temp_out_names= utility::split( base_name );
		utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
		base_name = out_name.base();
		InterfaceHit ih(cb_score, 0.0,0.0,0.0, 0 );
		std::string outname = base_name+option[ out::suffix ]()+"_"+right_string_of( 1, 8, '0' )+".pdb";
		dump_transformed_pdb( pose, ih, urs, outname, base_name );

		return;
	}

	int parallel_i = option[ run_i ]();
	int parallel_j = option[ run_j ]();
	for ( int ctr=(int)rot_lb; ctr<=(int)rot_ub ; ++ctr ) {
		if ( ctr%parallel_j != parallel_i%parallel_j ) continue;

		TR << "Rotation " << ctr << " of " << urs.nrots() << std::endl;
		urs.get(ctr, r_local);

		// interface_map:           stores the exact transformation and area of all interfaces > minintarea
		// ambiguous_interface_map: stores the index (in 'rts') of all interfaces with a sum > minintarea
		// p1_interface_map:        self interactions, independent of translation
		FArray3D< Real > sum_interface_area;
		//FArray3D< utility::vector1<Size> > ambiguous_interface_map;
		utility::vector1<SingleInterface> p1_interface_map;

		//ambiguous_interface_map.dimension( grid_[0] , grid_[1] , grid_[2] );
		sum_interface_area.dimension( grid_[0] , grid_[1] , grid_[2] ); sum_interface_area=0;
		Real self_ca = resample_maps_and_get_self( rho_ca, rho_cb, r_local, sg_, r_rho_ca, r_rho_cb, p1_interface_map );

		if ( debug_||debug_exact_ ) {
			std::ostringstream oss1; oss1 << "rot"<<ctr<<".mrc";
			writeMRC( r_rho_ca, oss1.str(), false, true );
			std::ostringstream oss2; oss2 << "rot"<<ctr<<".pdb";
			dump_transformed_pdb( pose, InterfaceHit( 0.,0.,0.,0., ctr ), urs, oss2.str(),"null_base_name" );
		}

		if ( self_ca >= maxclash ) {
			TR << "   self clashing!" << std::endl;
			if ( debug_ ) writeMRC( r_rho_ca, "clash.mrc" );
			continue; // next rotation
		}

		// P1 interactions are OK, now compute the rest
		FArray3D<Real> collision_map, ex_collision_map;
		collision_map.dimension( grid_[0] , grid_[1] , grid_[2] );
		collision_map=self_ca;

		Real sum_interface_p1 = 0;
		for ( int i=1; i<=(int)p1_interface_map.size(); ++i ) {
			sum_interface_p1 += p1_interface_map[i].cb_overlap_;
		}
		sum_interface_area = sum_interface_p1;

		if ( debug_exact_ ) {
			ex_collision_map.dimension( grid_[0] , grid_[1] , grid_[2] );
			ex_collision_map=self_ca;
		}

		// do the convolution with each symmop to find configurations that:
		//   (a) are not clashing
		//   (b) _might_ be fully connected in the lattice
		for ( int s=sym_lb; s<=(int)sym_ub; ++s ) {   // s==1 is always identity
			numeric::xyzMatrix<Real> s_i = rts[s].get_rotation();
			numeric::xyzMatrix<Real> s_inv = numeric::inverse(s_i);
			numeric::xyzVector<Real> t_inv = s_inv*(-rts[s].get_translation());

			transform_map( r_rho_ca, s_inv, t_inv, working_s);

			// all the magic is in here
			do_convolution( r_rho_ca, working_s, s_i, conv_out);

			for ( int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z ) {
				for ( int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y ) {
					for ( int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x ) {
						collision_map(x,y,z) += conv_out(x,y,z)*voxel_volume_;
					}
				}
			}

			// debug: exact CA fft
			if (  debug_exact_ ) {
				for ( int sz=(int)ccIndexLow[2]-1; sz<=(int)ccIndexHigh[2]-1; ++sz ) {
					for ( int sy=(int)ccIndexLow[1]-1; sy<=(int)ccIndexHigh[1]-1; ++sy ) {
						for ( int sx=(int)ccIndexLow[0]-1; sx<=(int)ccIndexHigh[0]-1; ++sx ) {
							ex_collision_map(sx+1,sy+1,sz+1) += get_clash_score_exact( numeric::xyzVector<int>(sx,sy,sz), rts[s].get_rotation(), rts[s].get_translation(), r_rho_ca);
						}
					}
				}
			}

			//int translation_count=1;
			// do CB fft
			transform_map( r_rho_cb, s_inv, t_inv, working_s);
			do_convolution( r_rho_cb, working_s, s_i, conv_out);
			for ( int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z ) {
				for ( int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y ) {
					for ( int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x ) {
						sum_interface_area(x,y,z) += conv_out(x,y,z)*voxel_volume_;
						translation_count=translation_count+1;
					}
				}
			}
		}

		int total_possible=translation_count*urs.nrots();
		TR<< "Size of search space: " << total_possible<<std::endl;

		if ( debug_ || debug_exact_ ) {
			std::ostringstream oss; oss << "collisionmap_"<<ctr<<".mrc";
			FArray3D<Real> collision_map_dump = collision_map;
			for ( int i=0; i<(int)Npoints; ++i ) collision_map_dump[i] = /*maxclash-*/collision_map[i];
			writeMRC( collision_map_dump, oss.str() );

			if ( debug_exact_ ) {
				std::ostringstream oss2; oss2 << "ex_collisionmap_"<<ctr<<".mrc";
				for ( int i=0; i<(int)Npoints; ++i ) collision_map_dump[i] = /*maxclash-*/ex_collision_map[i];
				writeMRC( collision_map_dump, oss2.str() );

				// get correl
				Real x2=0,y2=0,xy=0,x=0,y=0, Ncc=0;
				for ( int sz=(int)ccIndexLow[2]-1; sz<=(int)ccIndexHigh[2]-1; ++sz ) {
					for ( int sy=(int)ccIndexLow[1]-1; sy<=(int)ccIndexHigh[1]-1; ++sy ) {
						for ( int sx=(int)ccIndexLow[0]-1; sx<=(int)ccIndexHigh[0]-1; ++sx ) {
							x2 += collision_map(sx+1,sy+1,sz+1) * collision_map(sx+1,sy+1,sz+1);
							y2 += ex_collision_map(sx+1,sy+1,sz+1) * ex_collision_map(sx+1,sy+1,sz+1);
							xy += collision_map(sx+1,sy+1,sz+1) * ex_collision_map(sx+1,sy+1,sz+1);
							x  += collision_map(sx+1,sy+1,sz+1);
							y  += ex_collision_map(sx+1,sy+1,sz+1);
							Ncc++;
						}
					}
				}
				Real correl = 1, slope = 0;
				Real sx = (Ncc*x2-x*x);
				Real sy = (Ncc*y2-y*y);
				if ( sx*sy > 0 ) {
					correl = (Ncc*xy - x*y) / std::sqrt( (sx) * (sy) );
					slope = correl * sx/sy;
				}
				TR << "correl = " << correl << "   scale = " << slope << std::endl;

				collision_map = ex_collision_map;
			}
		}

		// finally add nonclashing interfaces to the DB
		Real mininterfacesum_filter = option[crystdock::mininterfacesum]();
		for ( int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z ) {
			for ( int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y ) {
				for ( int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x ) {
					if ( collision_map(x,y,z) < maxclash && sum_interface_area(x,y,z) > mininterfacesum_filter ) {
						nnonclashing++;

						// get_interface_score populates iinfo
						//    then computes the weakest connection necessary to construct the lattice
						utility::vector1<SingleInterface> iinfo;
						numeric::xyzVector<Real> xyz((Real)x-1,(Real)y-1,(Real)z-1);
						xyz = i2c_*xyz;

						get_interfaces_allatom( pose, rts, r_local, xyz, iinfo );
						Real score_xyz = get_interface_score (pose, iinfo, rts );

						if ( score_xyz >= mininterface_ ) {
							nconnected++;
							IDB.add_interface( InterfaceHit( score_xyz, xyz[0],xyz[1],xyz[2], ctr ) );
						}
					}
				}
			}
		}
		if ( IDB.size()>0 ) {
			TR << IDB.size() << " of " << nnonclashing << " nonclashing and " << nconnected << " connected "
				<< " configurations; min_score = " << IDB.top().score << std::endl;
		} else {
			TR << IDB.size() << " of " << nnonclashing << " nonclashing configurations" << std::endl;
		}
	}

	utility::vector1< InterfaceHit > ih_vec;
	int nhits = IDB.size();
	for ( int i=1; i<=nhits; ++i ) {
		InterfaceHit ih = IDB.pop();
		ih_vec.push_back( ih );
	}
	std::reverse(ih_vec.begin(),ih_vec.end());


	Real pose_radius = get_radius_of_pose( pose );
	Real cluster_cutoff=option[crystdock::cluster_cutoff];
	utility::vector1< InterfaceHit > ih_vec_clustered;
	for ( int i=1; i<=nhits; ++i ) {
		int nclust = ih_vec_clustered.size();
		bool found_match=false;
		for ( int j=1; j<=nclust && !found_match; ++j ) {
			Real dist_ij = get_transform_distance( ih_vec[i], ih_vec_clustered[j] , urs, pose_radius);
			found_match = (dist_ij <= cluster_cutoff);
		}
		if ( !found_match ) {
			ih_vec_clustered.push_back( ih_vec[i] );
		}
	}

	int nhits_after_cluster = ih_vec_clustered.size();
	TR << "Have " << nhits_after_cluster << " after clustering." << std::endl;
	int n_output=std::min(nhits_after_cluster, (int)nmodels_);
	for ( int i=1; i<=n_output; ++i ) {
		InterfaceHit ih = ih_vec_clustered[i];
		TR << i << ": " << ih.score << " " << ih.x << " "  << ih.y << " "  << ih.z << " " << ih.rot_index  << " " << std::endl;

		// Treat tags as file names so that we put the number before the extension.
		std::string base_name = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();
		utility::vector1< std::string > temp_out_names= utility::split( base_name );
		utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
		base_name = out_name.base();
		std::string outname = base_name+option[ out::suffix ]()+"_"+right_string_of( i, 8, '0' )+".pdb";
		dump_transformed_pdb( pose, ih, urs, outname,base_name );
	}
}

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( MoverOP( new CrystDock() ) );

	// force some options
	option[ out::nooutput ].value(true);
	option[ in::preserve_crystinfo ].value(true);

	// main loop
	protocols::jd2::JobDistributor::get_instance()->go( seq );

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	try {
		NEW_OPT(crystdock::spacegroup, "spacegroup (no spaces)", "P23");
		NEW_OPT(crystdock::A, "unit cell A", 50);
		NEW_OPT(crystdock::B, "unit cell B", 50);
		NEW_OPT(crystdock::C, "unit cell C", 50);
		NEW_OPT(crystdock::alpha, "unit cell alpha", 90);
		NEW_OPT(crystdock::beta, "unit cell beta", 90);
		NEW_OPT(crystdock::gamma, "unit cell gamma", 90);
		NEW_OPT(crystdock::maxclash, "max allowed clashscore", 20);
		NEW_OPT(crystdock::mininterface, "min allowed interface area", 250);
		NEW_OPT(crystdock::mininterfacesum, "min interface area sum", 4000);
		NEW_OPT(crystdock::trans_step, "translational stepsize (A)", 1);
		NEW_OPT(crystdock::rot_step, "rotational stepsize (degrees) ([debug] 0 searches input rotation only)", 10);
		NEW_OPT(crystdock::nmodels, "number of models to output", 1000);
		NEW_OPT(crystdock::ssonly, "limit interface calcs to secstruct only", false);
		NEW_OPT(crystdock::rotnum, "[debug] only run a single rotation", 0);
		NEW_OPT(crystdock::symnum, "[debug] only run a single symmop", 0);
		NEW_OPT(crystdock::debug, "[debug] dump intermediate info", false);
		NEW_OPT(crystdock::debug_exact, "[debug] debug mode with exact (non-FFT) calculations (slow!)", false);
		NEW_OPT(crystdock::eval_native, "[debug] evaluate input structure without docking", false);
		NEW_OPT(crystdock::randomize_orientation, "randomize orientation of input", false);
		NEW_OPT(crystdock::random_rotate, "random fixed angle rotation (eval native only)", 0);
		NEW_OPT(crystdock::n_clashdist, "n_clashdist", 1.40);
		NEW_OPT(crystdock::ca_clashdist, "ca_clashdist", 2.00);
		NEW_OPT(crystdock::c_clashdist, "c_clashdist", 2.00);
		NEW_OPT(crystdock::o_clashdist, "o_clashdist", 1.30);
		NEW_OPT(crystdock::cb_clashdist, "cb_clashdist", 1.50);
		NEW_OPT(crystdock::sigwidth, "sigwidth", 18.00);
		NEW_OPT(crystdock::interfacedist, "interfacedistance", 5.50);
		NEW_OPT(crystdock::interface_sigwidth, "interface_sigwidth", 1.00);
		NEW_OPT(crystdock::cluster_cutoff, "cluster_cutoff", 2.00);
		NEW_OPT(crystdock::expand_rounds, "expand_rounds", 12);
		NEW_OPT(crystdock::compact, "output transformation matrix only", false);
		NEW_OPT(run_i, "parallelize i of j", 0);
		NEW_OPT(run_j, "parallelize i of j", 1);

		devel::init( argc, argv );
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
