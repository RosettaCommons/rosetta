// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/UniformRotationSampler.cc
/// @author Frank DiMaio
// Uniformly sample rotation space
//   * use idea from A. Yershova, et al, International Journal of Robotics Research, 2009.
//       to get SO3 samples from S1+S2 uniform sampling
//   * use icosahedral embedding to generate approximately uniform coverage of S2


#ifndef INCLUDED_numeric_uniformrotationsampler_hh
#define INCLUDED_numeric_uniformrotationsampler_hh

#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/Quaternion.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/types.hh>

namespace numeric {

inline double urs_norm4(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

inline platform::Real
urs_R2ang(numeric::xyzMatrix<Real> R) {
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
	Real r = urs_norm4(q0, q1, q2, q3);
	q0 /= r; q1 /= r; q2 /= r; q3 /= r;
	Real angle = 2*conversions::degrees( std::fabs(acos(q0)) );

	return (angle);
}


struct urs_Quat {
	Real x_,y_,z_,w_;
	urs_Quat() {
		x_=0; y_=0; z_=0; w_=1;
	}
	urs_Quat( Real x_in,Real y_in,Real z_in,Real w_in) {
		x_=x_in; y_=y_in; z_=z_in; w_=w_in;
	}

	// construct from a rotation matrix
	// should be stable
	urs_Quat( numeric::xyzMatrix<Real> const& R) {
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
			return numeric::xyzMatrix<platform::Real>::rows(1,0,0, 0,1,0, 0,0,1);
		}

		Real xx = x_*x_, xy = x_*y_, xz = x_*z_, xw = x_*w_;
		Real yy = y_*y_, yz = y_*z_, yw = y_*w_;
		Real zz = z_*z_, zw = z_*w_;
		//Real ww = w_*w_;

		return numeric::xyzMatrix<platform::Real>::rows(
			1 - 2 * ( yy+zz ) ,     2 * ( xy-zw ) ,     2 * ( xz+yw ) ,
			2 * ( xy+zw ) , 1 - 2 * ( xx+zz ) ,     2 * ( yz-xw ) ,
			2 * ( xz-yw ) ,     2 * ( yz+xw ) , 1 - 2 * ( xx+yy ) );
	}
}; // urs_Quat


class UniformRotationSampler {
private:
	utility::vector1< urs_Quat > rotlist_;
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
		utility::vector1<platform::Real> cos_psi_over_2(S1subs), sin_psi_over_2(S1subs);
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
				rotlist_[counter++] = urs_Quat(x,y,z,w);
			}
		}
	}


	void
	remove_redundant( numeric::xyzMatrix<numeric::Real> R0 ) {
		utility::vector1<bool> tokeep(rotlist_.size(), true);
		utility::vector1<numeric::xyzMatrix<Real> > Rs(rotlist_.size());
		utility::vector1<numeric::xyzMatrix<Real> > Rinvs(rotlist_.size());

		for ( int i=1; i<=(int)rotlist_.size(); ++i ) {
			Rs[i] = rotlist_[i].asR();
			Rinvs[i] = numeric::inverse( Rs[i] );
		}


		for ( int i=1; i<=(int)rotlist_.size(); ++i ) {
			if ( !tokeep[i] ) { continue; }
			numeric::xyzMatrix<Real> SR = R0*Rs[i];
			for ( int j=i+1; j<=(int)rotlist_.size(); ++j ) {
				if ( !tokeep[j] ) { continue; }

				Real ang_ij = urs_R2ang( SR*Rinvs[j] );
				if ( ang_ij<0.65*theta_ ) {
					tokeep[j]=false;
				}
			}
		}

		utility::vector1< urs_Quat > rotlist_old = rotlist_;
		rotlist_.clear();
		for ( int i=1; i<=(int)rotlist_old.size(); ++i ) {
			if ( tokeep[i] ) rotlist_.push_back(rotlist_old[i]);
		}

		std::cout << "Trimmed rot list has " << rotlist_.size() << " rotations (was " << rotlist_old.size() << ")\n";
	}
}; // class UniformRotationSampler

} // namespace numeric

#endif
