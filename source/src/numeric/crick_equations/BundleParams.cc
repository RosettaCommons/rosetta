// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file BundleParams.cc.
/// @brief Functions implementing the Crick equations for a helical bundle (a helix of helices).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//C++ headers
#include <math.h>

// Unit headers
#include <numeric/types.hh>
#include <numeric/constants.hh>

//Project headers
#include <numeric/crick_equations/BundleParams.hh>

using namespace std;

namespace numeric {
namespace crick_equations {

/// @brief Calculates P, the repeat distance of the major helix in the z-direction.
/// @details Sets "failed" to TRUE if and only if the value of P could not be calculated.
Real PP (
	Real const &r0,
	Real const &omega0,
	Real const &z1,
	bool &failed
) {
	Real const val = pow(z1/omega0,2) - pow(r0,2) ;
	if ( val < 0 ) {
		failed = true;
		return 0.0;
	}
	failed=false;
	return numeric::constants::d::pi_2*sqrt( val );
}

/// @brief Calculates P*omega0 (the repeat distance of the major helix in the z-direction times the major helix
/// turn per residue).
/// @details Sets "failed" to TRUE if and only if the value of P*omega0 could not be calculated.
Real P_omega0 (
	Real const &r0,
	Real const &omega0,
	Real const &z1,
	bool &failed
) {
	Real const val = pow(z1,2) - pow(r0*omega0,2) ;
	if ( val < 0 ) {
		failed = true;
		return 0.0;
	}
	failed=false;
	return numeric::constants::d::pi_2*sqrt( val );
}

/// @brief Calculates alpha, the tilt angle of the minor helix.
///
Real ALPHA (
	Real const &r0,
	Real const &omega0,
	Real const &Pomega0 //P times omega0
) {
	//return atan( numeric::constants::d::pi_2*r0/P );
	return atan2( numeric::constants::d::pi_2*r0*omega0, Pomega0 );
}

/// @brief Calculates cos(omega*t+delta_omega)
///
Real COSFXN( Real const &t, Real const &omega, Real const&delta_omega ) {
	return cos(omega*t+delta_omega);
}

/// @brief Calculates the norm of the gradient vector
/// @details Needed for adding the small delta_z1 offsets, which are in the direction
/// of (dx/dt, dy/dt, dz/dt) of the major helix.
Real gradnorm(
	Real const &r0,
	Real const &omega0,
	Real const &Pomega0, //P times omega0
	Real const &s0,
	Real const &c0
) {
	return sqrt( pow(-r0*omega0*s0,2) + pow(r0*omega0*c0,2) + pow(Pomega0/numeric::constants::d::pi_2,2) );
}

/// @brief Calculates sin(omega*t+delta_omega)
///
Real SINFXN( Real const &t, Real const &omega, Real const&delta_omega ) {
	return sin(omega*t+delta_omega);
}

/// @brief Calculates the x-position on the helix of helices
/// given the Crick parameters for the bundle.
/// @details Returns failed=true if calculation fails, failed=false otherwise.
Real X_BUNDLE (
	Real const &t,
	Real const &r0,
	Real const &omega0,
	Real const &delta_omega0,
	//Real const &P,
	Real const &r1,
	Real const &omega1,
	Real const &z1,
	Real const &delta_omega1,
	Real const &delta_z1,
	bool &failed
) {
	//Real const P = PP(r0,omega0,z1, failed);
	Real const Pomega0 = P_omega0(r0,omega0,z1,failed);
	if ( failed ) return 0.0;
	Real const alpha=ALPHA(r0,omega0,Pomega0);
	Real const c0 = COSFXN( t, omega0, delta_omega0 );
	Real const c1 = COSFXN( t, omega1, delta_omega1 );
	Real const s0 = SINFXN( t, omega0, delta_omega0 );
	Real const s1 = SINFXN( t, omega1, delta_omega1 );
	return r0*c0+r1*c0*c1-r1*cos(alpha)*s0*s1 - (r0*omega0*s0/gradnorm(r0,omega0,Pomega0,s0,c0))*delta_z1;
}

/// @brief Calculates the y-position on the helix of helices
/// given the Crick parameters for the bundle.
/// @details Returns failed=true if calculation fails, failed=false otherwise.
Real Y_BUNDLE (
	Real const &t,
	Real const &r0,
	Real const &omega0,
	Real const &delta_omega0,
	//Real const &P,
	Real const &r1,
	Real const &omega1,
	Real const &z1,
	Real const &delta_omega1,
	Real const &delta_z1,
	bool &failed
) {
	//Real const P = PP(r0,omega0,z1,failed);
	Real const Pomega0 = P_omega0(r0,omega0,z1,failed);
	if ( failed ) return 0.0;
	Real const alpha=ALPHA(r0,omega0,Pomega0);
	Real const c0 = COSFXN( t, omega0, delta_omega0 );
	Real const c1 = COSFXN( t, omega1, delta_omega1 );
	Real const s0 = SINFXN( t, omega0, delta_omega0 );
	Real const s1 = SINFXN( t, omega1, delta_omega1 );
	return r0*s0+r1*s0*c1+r1*cos(alpha)*c0*s1 + (r0*omega0*c0/gradnorm(r0,omega0,Pomega0,s0,c0))*delta_z1;
}

/// @brief Calculates the z-position on the helix of helices
/// given the Crick parameters for the bundle.
/// @details Returns failed=true if calculation fails, failed=false otherwise.
Real Z_BUNDLE (
	Real const &t,
	Real const &r0,
	Real const &omega0,
	Real const &delta_omega0,
	//Real const &P,
	Real const &r1,
	Real const &omega1,
	Real const &z1,
	Real const &delta_omega1,
	Real const &delta_z1,
	bool &failed
) {
	//Real const P = PP(r0,omega0,z1,failed);
	Real const Pomega0 = P_omega0(r0,omega0,z1,failed);
	if ( failed ) return 0.0;
	Real const alpha=ALPHA(r0,omega0,Pomega0);
	Real const c0 = COSFXN( t, omega0, delta_omega0 );
	//Real const c1 = COSFXN( t, omega1, delta_omega1 );
	Real const s0 = SINFXN( t, omega0, delta_omega0 );
	Real const s1 = SINFXN( t, omega1, delta_omega1 );
	return (Pomega0*t/numeric::constants::d::pi_2-r1*sin(alpha)*s1 + (Pomega0/numeric::constants::d::pi_2/gradnorm(r0,omega0,Pomega0,s0,c0))*delta_z1);
}

/// @brief Calculate the x,y,z coordinates of a point on the helix of
/// helices given the Crick parameters for the bundle.
/// @details Not quite as efficient as it could be, but it probably
/// doesn't matter.  (If I really wanted to optimize this, I'd make sure
/// that c0, c1, s0, s1, and alpha were all calculated once rather than
/// thrice.)  Returns failed=true if calculation fails, failed=false otherwise.
xyzVector <Real> XYZ_BUNDLE (
	Real const &t,
	Real const &r0,
	Real const &omega0,
	Real const &delta_omega0,
	//Real const &P,
	Real const &r1,
	Real const &omega1,
	Real const &z1,
	Real const &delta_omega1,
	Real const &delta_z1,
	bool &failed
) {
	bool failedx=false, failedy=false, failedz=false;
	xyzVector < Real > returnvector(0,0,0);
	returnvector.x( X_BUNDLE( t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1, failedx ) );
	returnvector.y( Y_BUNDLE( t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1, failedy ) );
	returnvector.z( Z_BUNDLE( t, r0, omega0, delta_omega0, r1, omega1, z1, delta_omega1, delta_z1, failedz ) );
	failed=failedx||failedy||failedz;
	if ( failed ) returnvector.assign(0,0,0);
	return returnvector;
}

} //namespace crick_equations
} //namespace numeric
