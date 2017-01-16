// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_fiber_diffraction_xray_scattering_hh
#define INCLUDED_core_scoring_fiber_diffraction_xray_scattering_hh

#include <core/types.hh>

#include <cmath>
#include <string>

#ifdef WIN32
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#endif


// C++ Headers

namespace core {
namespace scoring {
namespace fiber_diffraction {


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
		return (c_ + a_[0]*exp(-b_[0]*S2/4) + a_[1]*exp(-b_[1]*S2/4) + a_[2]*exp(-b_[2]*S2/4) + a_[3]*exp(-b_[3]*S2/4) );
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
	OneGaussianScattering(int w, float s) {
		sigma_ = s;
		weight_ = w;
	}

	// rho = C*exp(-k*X^2)
	inline core::Real k( core::Real B, core::Real min_grid ) const {
		// s <= real-space sigma
		//    sig = sig+B/4;    % conv by B
		//    sig = max( sig, (reso).^2 );  % res limit
		//    D2 = (x-atom).^2;
		//    rho_c_alt2(ele,:) = grid * ele * sqrt(pi/sig) * exp(-(pi^2/sig).*D2);
		core::Real s = std::max( sigma_ + B/4 , 4*(min_grid*min_grid) ) ;
		return ( M_PI*M_PI/s );
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
	inline float s( ) const {
		return ( sigma_ );
	}

private:
	float sigma_;
	int weight_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


////////////
// Precomputed scattering
// fpd -- doesn't seem to make calculation much faster
//
// class AtomScattering {
// public:
//  AtomScattering();
//
//  // mask scattering
//  AtomScattering( core::Real mask, core::Real reso );
//
//  // atom scattering (single gaussian)
//  AtomScattering( core::Real a, core::Real B, core::Real mask, core::Real reso );
//
//  // TODO: atom scattering (crystallographic form-factors)
//  //AtomScattering( KromerMann f0, core::Real mask, core::Real reso );
//
//  // interp at a Cartesian point X
//  core::Real interp_linear( numeric::xyzVector< core::Real > const & X ) const;
//
// private:
//  // common initialization
//  void init( core::Real mask, core::Real reso );
//
//  ObjexxFCL::FArray3D< double > data;
//  numeric::xyzVector< core::Real > c2i, i2c;    // always construct orthogonal
//  numeric::xyzVector< core::Size > grid;
//  core::Real reso;
// };


// weight from scattering factors
OneGaussianScattering get_A( std::string elt );

// weight from scattering factors
KromerMann get_km( std::string elt );

// precomputed scattering
//const AtomScattering & get_scattering( std::string elt );

bool factorsLTE5(int X);
bool factorsLTE19(int X);
int findSampling(double MINSMP, int NMUL);
int findSampling5(double MINSMP, int NMUL);

}
}
}

#endif
