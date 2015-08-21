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
/// @author

#ifndef INCLUDED_core_scoring_electron_density_xray_scattering_hh
#define INCLUDED_core_scoring_electron_density_xray_scattering_hh

#include <core/types.hh>

#include <cmath>
#include <string>


#ifdef WIN32

#define _USE_MATH_DEFINES

#include <math.h>

#endif


// C++ Headers

namespace core {
namespace scoring {
namespace electron_density {


////////////
// 4-gaussian approximation of scattering coefficients
class KromerMann {
public:
	KromerMann() {
		a_[0] = a_[1] = a_[2] = a_[3] = 0;
		b_[0] = b_[1] = b_[2] = b_[3] = 0;
		c_ = 0;
	}
	KromerMann(float c, float a1, float a2, float a3, float a4, float b1, float b2, float b3, float b4) {
		a_[0] = a1; a_[1] = a2; a_[2] = a3; a_[3] = a4;
		b_[0] = b1; b_[1] = b2; b_[2] = b3; b_[3] = b4;
		c_ = c;
	}

	// scattering at reciprocal space distance^2
	inline core::Real f0( core::Real S2 ) {
		return (c_ + a_[0]*exp(b_[0]*S2/4) + a_[1]*exp(b_[1]*S2/4) + a_[2]*exp(b_[2]*S2/4) + a_[3]*exp(b_[3]*S2/4) );
	}

private:
	float a_[4];
	float b_[4];
	float c_;
};

////////////
// 1-gaussian real-space approximate scattering
class OneGaussianScattering {
public:
	OneGaussianScattering() {
		weight_ = 0;
		sigma_ = 3.0;
	}
	OneGaussianScattering(float w, float s) {
		sigma_ = s;
		weight_ = w;
	}

	inline core::Real B( core::Real k ) const {
		core::Real s = M_PI*M_PI/k;
		core::Real sigma_eff = sigma_;
		core::Real B = ( 4*( s - sigma_eff ) );

		// smooth to flat at B==0
		if ( B < 0 ) B = 0;
		else if ( B<10 ) B = sqrt(10*B);
		return ( B );
	}

	inline core::Real k( core::Real B , core::Real lim=600) const {
		core::Real sigma_eff = sigma_;
		core::Real B_eff = B;
		if ( B<0 ) B_eff = 0;
		else if ( B<1 ) B_eff = B*B;
		else if ( B>lim-10.0 ) B_eff = lim-(1.0/10.0)*(lim-B)*(lim-B);
		else if ( B>lim ) B_eff = lim;
		core::Real s = sigma_eff + B_eff/4;
		core::Real k = M_PI*M_PI/s;
		return k;
	}

	// calculate dK/dB at a given resolution
	inline core::Real dk( core::Real B , core::Real lim=600) const {
		core::Real sigma_eff = sigma_;

		core::Real B_eff = B;
		if ( B<0 || B>lim ) return 0;
		else if ( B<1 ) B_eff = B*B;
		else if ( B>lim-10.0 ) B_eff = lim-(1.0/10.0)*(lim-B)*(lim-B);

		core::Real s = sigma_eff + B_eff/4;
		core::Real dkdb = -M_PI*M_PI/(4*s*s);

		if ( B<1 ) dkdb *= 2*B;
		if ( B>lim-10.0 ) dkdb *= (1.0/5.0)*(lim-B);

		return dkdb;
	}

	// rho = C*exp(-k*X^2)
	// get scale factor, given k
	inline core::Real C( core::Real k ) const {
		core::Real C = pow(k, 1.5);
		return ( C*weight_ );
	}

	inline int a( ) const {
		return ( weight_ );
	}

private:
	float sigma_;
	float weight_;
};


// weight from scattering factors
OneGaussianScattering get_A( std::string elt );

// weight from scattering factors
KromerMann get_km( std::string elt );

bool factorsLTE5(int X);
bool factorsLTE19(int X);
int findSampling(double MINSMP, int NMUL);
int findSampling5(double MINSMP, int NMUL);

}
}
}

#endif
