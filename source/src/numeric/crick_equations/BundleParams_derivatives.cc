// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/crick_equations/BundleParams_derivatives.cc
/// @brief  Analytical derivatives of the Crick equations with respect to bundle parameters.
/// @details For the epsilon=1 (circular cross-section) case, the Crick equations are:
///
///   Pω0 = 2π · √(z1² - (r0·ω0)²)
///   α   = atan2(2π·r0·ω0, Pω0)
///   gn  = z1  (simplifies when epsilon=1)
///
///   X = r0·c0 + r1·c1·c0 - r1·cos(α)·s0·s1 - (r0·ω0·s0/gn)·δz1
///   Y = r0·s0 + r1·c1·s0 + r1·cos(α)·c0·s1 + (r0·ω0·c0/gn)·δz1
///   Z = (Pω0·t/2π) - r1·sin(α)·s1 + (Pω0/2π/gn)·δz1
///
/// where c0=cos(ω0·t+δω0), s0=sin(ω0·t+δω0), c1=cos(ω1·t+δω1), s1=sin(ω1·t+δω1).
///
/// @author Andy Watkins

#include <numeric/crick_equations/BundleParams_derivatives.hh>
#include <numeric/constants.hh>

#include <cmath>

namespace numeric {
namespace crick_equations {

CrickDerivatives compute_crick_derivatives(
	Real const t,
	Real const r0,
	Real const omega0,
	Real const delta_omega0,
	Real const r1,
	Real const omega1,
	Real const z1,
	Real const delta_omega1,
	Real const delta_z1,
	bool & failed
) {
	CrickDerivatives derivs;
	failed = false;

	Real const twopi = numeric::constants::d::pi_2;

	// Intermediate quantities
	Real const r0w0 = r0 * omega0;
	Real const r0w0_sq = r0w0 * r0w0;
	Real const z1_sq = z1 * z1;
	Real const val = z1_sq - r0w0_sq;
	if ( val < 0 ) {
		failed = true;
		return derivs;
	}
	Real const sqrt_val = std::sqrt( val );
	Real const Pw0 = twopi * sqrt_val;              // P*omega0
	Real const Pw0_over_twopi = sqrt_val;            // Pw0 / 2pi = sqrt(z1^2 - (r0*w0)^2)
	Real const gn = z1;                              // gradnorm = z1 when epsilon=1

	// alpha and its trig
	Real const alpha = std::atan2( twopi * r0w0, Pw0 );  // = atan2(r0*w0, sqrt_val)
	Real const cos_alpha = std::cos( alpha );
	Real const sin_alpha = std::sin( alpha );

	// Trig basis
	Real const c0 = std::cos( omega0 * t + delta_omega0 );
	Real const s0 = std::sin( omega0 * t + delta_omega0 );
	Real const c1 = std::cos( omega1 * t + delta_omega1 );
	Real const s1 = std::sin( omega1 * t + delta_omega1 );

	// Ratio used in delta_z1 terms
	Real const w0_over_gn = omega0 / gn;           // omega0 / z1
	Real const r0w0_over_gn = r0w0 / gn;           // r0*omega0 / z1
	Real const Pw0_over_twopi_gn = Pw0_over_twopi / gn;  // sqrt(z1^2-(r0*w0)^2) / z1

	// ========================
	// dXYZ / d(delta_omega0)
	// ========================
	// delta_omega0 only appears in c0, s0: dc0/ddw0 = -s0, ds0/ddw0 = c0
	{
		Real const dX = -r0 * s0 - r1 * c1 * s0 - r1 * cos_alpha * c0 * s1 - r0w0_over_gn * c0 * delta_z1;
		Real const dY =  r0 * c0 + r1 * c1 * c0 - r1 * cos_alpha * s0 * s1 - r0w0_over_gn * s0 * delta_z1;
		Real const dZ = 0.0;
		derivs.dXYZ_ddelta_omega0 = xyzVector<Real>( dX, dY, dZ );
	}

	// ========================
	// dXYZ / d(delta_omega1)
	// ========================
	// delta_omega1 only appears in c1, s1: dc1/ddw1 = -s1, ds1/ddw1 = c1
	{
		Real const dX = -r1 * s1 * c0 - r1 * cos_alpha * s0 * c1;
		Real const dY = -r1 * s1 * s0 + r1 * cos_alpha * c0 * c1;
		Real const dZ = -r1 * sin_alpha * c1;
		derivs.dXYZ_ddelta_omega1 = xyzVector<Real>( dX, dY, dZ );
	}

	// ========================
	// dXYZ / d(delta_t)
	// ========================
	// delta_t shifts t → t + delta_t. So dXYZ/d(delta_t) = dXYZ/dt.
	// This is the tangent vector to the helix at position t.
	{
		// dc0/dt = -omega0*s0, ds0/dt = omega0*c0
		// dc1/dt = -omega1*s1, ds1/dt = omega1*c1
		Real const dX = -r0 * omega0 * s0
			+ r1 * (-omega1 * s1 * c0 - omega0 * s0 * c1)
			- r1 * cos_alpha * (omega0 * c0 * s1 + omega1 * s0 * c1)
			- r0w0_over_gn * omega0 * c0 * delta_z1;
		Real const dY = r0 * omega0 * c0
			+ r1 * (-omega1 * s1 * s0 + omega0 * c0 * c1)
			+ r1 * cos_alpha * (-omega0 * s0 * s1 + omega1 * c0 * c1)
			- r0w0_over_gn * omega0 * s0 * delta_z1;
		Real const dZ = Pw0_over_twopi - r1 * sin_alpha * omega1 * c1;
		derivs.dXYZ_ddelta_t = xyzVector<Real>( dX, dY, dZ );
	}

	// ========================
	// dXYZ / d(r0)
	// ========================
	// r0 appears directly AND through Pw0 and alpha.
	// dPw0/dr0 = 2pi * d/dr0[sqrt(z1^2-(r0*w0)^2)] = 2pi * (-r0*w0^2) / sqrt_val = -twopi*r0*w0^2/sqrt_val
	// d(alpha)/dr0: alpha = atan2(twopi*r0*w0, Pw0)
	//   Let u = twopi*r0*w0, v = Pw0 = twopi*sqrt_val
	//   d(atan2(u,v))/dr0 = (v*du/dr0 - u*dv/dr0) / (u^2+v^2)
	//   du/dr0 = twopi*w0
	//   dv/dr0 = twopi*(-r0*w0^2)/sqrt_val = dPw0/dr0
	//   u^2+v^2 = (twopi*r0*w0)^2 + (twopi*sqrt_val)^2 = twopi^2 * z1^2
	//   So d(alpha)/dr0 = [Pw0*twopi*w0 - twopi*r0*w0*(-twopi*r0*w0^2/sqrt_val)] / (twopi^2*z1^2)
	//     = twopi*w0*[Pw0 + twopi*r0^2*w0^2/sqrt_val] / (twopi^2*z1^2)
	//     = w0*[sqrt_val + r0^2*w0^2/sqrt_val] / (twopi*z1^2)
	//     = w0*z1^2/(sqrt_val*twopi*z1^2)
	//     = w0/(twopi*sqrt_val)
	{
		Real const dPw0_dr0 = -twopi * r0 * omega0 * omega0 / sqrt_val;
		Real const dalpha_dr0 = omega0 / ( twopi * sqrt_val );
		Real const dcos_alpha_dr0 = -sin_alpha * dalpha_dr0;
		Real const dsin_alpha_dr0 =  cos_alpha * dalpha_dr0;

		// d(r0*c0)/dr0 = c0
		// d(r0*w0*s0/gn)/dr0 = w0*s0/gn  (gn=z1 independent of r0)
		Real const dX = c0 + r1 * 0.0 - r1 * dcos_alpha_dr0 * s0 * s1 - w0_over_gn * s0 * delta_z1;
		Real const dY = s0 + r1 * 0.0 + r1 * dcos_alpha_dr0 * c0 * s1 + w0_over_gn * c0 * delta_z1;
		Real const dZ = (dPw0_dr0 * t / twopi) - r1 * dsin_alpha_dr0 * s1 + (dPw0_dr0 / twopi / gn) * delta_z1;
		derivs.dXYZ_dr0 = xyzVector<Real>( dX, dY, dZ );
	}

	// ========================
	// dXYZ / d(omega0)
	// ========================
	// omega0 appears through c0, s0 (via t*omega0+delta_omega0), Pw0, alpha, and the delta_z1 terms.
	// dc0/domega0 = -t*s0, ds0/domega0 = t*c0
	// dPw0/domega0 = 2pi * d/domega0[sqrt(z1^2-r0^2*w0^2)] = 2pi*(-r0^2*w0)/sqrt_val
	// d(alpha)/domega0 = r0/(twopi*sqrt_val)  [analogous derivation to dr0 case]
	{
		Real const dPw0_dw0 = -twopi * r0 * r0 * omega0 / sqrt_val;
		Real const dalpha_dw0 = r0 / ( twopi * sqrt_val );
		Real const dcos_alpha_dw0 = -sin_alpha * dalpha_dw0;
		Real const dsin_alpha_dw0 =  cos_alpha * dalpha_dw0;

		// Direct trig derivatives from omega0 appearing in the argument of c0, s0:
		Real const dc0 = -t * s0;
		Real const ds0 =  t * c0;

		// X = r0*c0 + r1*c1*c0 - r1*cos(alpha)*s0*s1 - (r0*w0*s0/gn)*dz1
		Real const dX = r0 * dc0
			+ r1 * c1 * dc0
			- r1 * (dcos_alpha_dw0 * s0 + cos_alpha * ds0) * s1
			- (r0 * s0 / gn + r0w0_over_gn * ds0) * delta_z1;

		Real const dY = r0 * ds0
			+ r1 * c1 * ds0
			+ r1 * (dcos_alpha_dw0 * c0 + cos_alpha * dc0) * s1
			+ (r0 * c0 / gn + r0w0_over_gn * dc0) * delta_z1;

		Real const dZ = (dPw0_dw0 * t / twopi)
			- r1 * dsin_alpha_dw0 * s1
			+ (dPw0_dw0 / twopi / gn) * delta_z1;

		derivs.dXYZ_domega0 = xyzVector<Real>( dX, dY, dZ );
	}

	return derivs;
}

} // namespace crick_equations
} // namespace numeric
