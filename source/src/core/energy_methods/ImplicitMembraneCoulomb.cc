// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/ImplicitMembraneCoulomb.cc
/// @brief Minimal class for computing the depth- and membrane-dependent electrostatics energy
/// @author rfalford12 (rfalford12@gmail.com)

// Project headers:
#include <core/energy_methods/ImplicitMembraneCoulomb.hh>

// Basic headers:
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/cubic_polynomial.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "core.energy_methods.ImplicitMembraneCoulomb" );

namespace core {
namespace energy_methods {

/// @brief Default constructor.
ImplicitMembraneCoulomb::ImplicitMembraneCoulomb():
	utility::VirtualBase(),
	max_dis_( 5.5 ),
	max_dis2_( 0.0 ),
	C0_( 0.0 ),
	C1_( 0.0 )
{
	initialize();
}

/// @brief Copy constructor.  Keep default unless deep copying is needed (and in that case,
/// consider using DeepCopyOPs.)
ImplicitMembraneCoulomb::ImplicitMembraneCoulomb( ImplicitMembraneCoulomb const & )=default;

/// @brief Destructor.
ImplicitMembraneCoulomb::~ImplicitMembraneCoulomb(){}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
ImplicitMembraneCoulombOP
ImplicitMembraneCoulomb::clone() const {
	return utility::pointer::make_shared< ImplicitMembraneCoulomb >( *this );
}

/// @brief Initialize max distance and coulomb constants
void
ImplicitMembraneCoulomb::initialize() {

	// set min & max distance
	max_dis_ = 5.5;
	max_dis2_ = max_dis_ * max_dis_;
	min_dis_ = 1.45;
	min_dis2_ = min_dis_ * min_dis_;

	// Set Constants
	C0_ = 322.0637;
	C1_ = C0_;

	dEfac_ = -1.0 * C0_;


	// Setup polynomials to fix derivative discontinuities
	low_poly_start_ = min_dis_ - 0.25;
	low_poly_end_ = min_dis_ + 0.25;
	low_poly_start2_ = low_poly_start_ * low_poly_start_;
	low_poly_end2_ = low_poly_end_ * low_poly_end_;

	hi_poly_start_ = max_dis_ - 1.0;
	hi_poly_end_ = max_dis_;
	hi_poly_start2_ = hi_poly_start_ * hi_poly_start_;
	hi_poly_end2_ = hi_poly_end_ * hi_poly_end_;


}

numeric::CubicPolynomial
ImplicitMembraneCoulomb::compute_lowpoly(core::Real const fi, core::Real const fj) const{       // scope low polynomial
	using namespace numeric::interpolation::spline;
	Real low_poly_end_score(0.0), low_poly_end_deriv(0.0);


	core::Real C2(0.0), min_dis_score(0.0);
	C2 = C1_ / (max_dis_ * compute_depth_and_bilayer_dep_dielectric( fi, fj, max_dis_ ));
	min_dis_score = compute_min_dis_score(fi, fj);

	//eps_low is dependent on fij too
	Real eps_low = compute_depth_and_bilayer_dep_dielectric( fi, fj, low_poly_end_ );
	Real deps_low = compute_deps_dr( fi, fj, low_poly_end_ );
	low_poly_end_score = C1_ / (low_poly_end_*eps_low) - C2;
	low_poly_end_deriv = -(C0_*(eps_low + low_poly_end_*deps_low))/(low_poly_end_*low_poly_end_*eps_low*eps_low);

	SplineGenerator gen_low_poly(
		low_poly_start_, min_dis_score, 0,
		low_poly_end_, low_poly_end_score, low_poly_end_deriv );
	InterpolatorOP interp_low( gen_low_poly.get_interpolator() );
	SimpleInterpolatorOP sinterp_low = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_low );
	if ( ! sinterp_low ) {
		utility_exit_with_message( "Hack Elec created non-simple-interpolator in initialize()" );
	}
	numeric::SplineParameters low_sp;
	low_sp.ylo  = sinterp_low->y()[ 1 ];
	low_sp.yhi  = sinterp_low->y()[ 2 ];
	low_sp.y2lo = sinterp_low->ddy()[ 1 ];
	low_sp.y2hi = sinterp_low->ddy()[ 2 ];
	numeric::CubicPolynomial low_poly(numeric::cubic_polynomial_from_spline( low_poly_start_, low_poly_end_, low_sp ));

	return(low_poly);
}

core::Real
ImplicitMembraneCoulomb::compute_min_dis_score(core::Real const fi,core::Real const fj ) const{
	core::Real C2(0.0);
	C2 = C1_ / (max_dis_ * compute_depth_and_bilayer_dep_dielectric( fi, fj, max_dis_ ));

	return( C1_ / (min_dis_ * compute_depth_and_bilayer_dep_dielectric( fi, fj, min_dis_ )) - C2);

}

numeric::CubicPolynomial
ImplicitMembraneCoulomb::compute_hipoly(core::Real const fi, core::Real const fj) const{
	// scope hi polynomial
	using namespace numeric::interpolation::spline;
	Real hi_poly_start_score(0.0), hi_poly_start_deriv(0.0);

	core::Real C2(0.0);
	C2 = C1_ / (max_dis_ * compute_depth_and_bilayer_dep_dielectric( fi, fj, max_dis_ ));


	Real eps_hi = compute_depth_and_bilayer_dep_dielectric( fi, fj, hi_poly_start_ );
	Real deps_hi = compute_deps_dr( fi, fj, hi_poly_start_ );
	hi_poly_start_score = C1_ /(hi_poly_start_*eps_hi) - C2;
	hi_poly_start_deriv = -(C0_*(eps_hi + hi_poly_start_*deps_hi))/(hi_poly_start_*hi_poly_start_*eps_hi*eps_hi);

	SplineGenerator gen_hi_poly(
		hi_poly_start_, hi_poly_start_score, hi_poly_start_deriv,
		hi_poly_end_, 0, 0 );

	InterpolatorOP interp_hi( gen_hi_poly.get_interpolator() );
	SimpleInterpolatorOP sinterp_hi = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_hi );
	if ( ! sinterp_hi ) {
		utility_exit_with_message( "Hack Elec created non-simple-interpolator in initialize()" );
	}
	numeric::SplineParameters hi_sp;
	hi_sp.ylo  = sinterp_hi->y()[ 1 ];
	hi_sp.yhi  = sinterp_hi->y()[ 2 ];
	hi_sp.y2lo = sinterp_hi->ddy()[ 1 ];
	hi_sp.y2hi = sinterp_hi->ddy()[ 2 ];

	numeric::CubicPolynomial hi_poly( numeric::cubic_polynomial_from_spline( hi_poly_start_, hi_poly_end_, hi_sp ) );
	return(hi_poly);
}


/// @brief Compute depth- and bilayer-dependent dielectric constant
core::Real
ImplicitMembraneCoulomb::compute_depth_and_bilayer_dep_dielectric(
	core::Real const fi, /* fractional hydration of atom i */
	core::Real const fj, /* fractional hydration of atom j */
	core::Real const d /* Distance between atoms i and j */
) const {

	core::Real f_ij = sqrt(fi*fj);
	core::Real D0_sol( 6.0 );
	core::Real D_sol( 80.0 );

	core::Real D0_mp( 3.0 );
	core::Real D_mp( 10.0 );
	core::Real slope( 0.40 );
	core::Real eps_sol(sigmoidal_eps(D0_sol, D_sol, slope, d));
	core::Real eps_mp(sigmoidal_eps(D0_mp, D_mp, slope, d));

	core::Real epsilon = ((1.0 - f_ij)*eps_mp) + (f_ij*eps_sol);
	return epsilon;
}

core::Real
ImplicitMembraneCoulomb::sigmoidal_eps( core::Real const sigmoid_D0,/*upper limit of sigmoid*/
	core::Real const sigmoid_D, /*lower limit of sigmoid*/
	core::Real const sigmoid_s,/*slope of sigmoid*/
	core::Real const d) const{

	return (sigmoid_D
		- 0.5 * (sigmoid_D - sigmoid_D0)
		* (2.0 + 2.0*d*sigmoid_s + d*d*sigmoid_s*sigmoid_s)
		* std::exp (-d*sigmoid_s) );
}

core::Real
ImplicitMembraneCoulomb::sigmoidal_deps_dr( core::Real const sigmoid_D0,/*upper limit of sigmoid*/
	core::Real const sigmoid_D, /*lower limit of sigmoid*/
	core::Real const sigmoid_s,/*slope of sigmoid*/
	core::Real const d) const{

	return ( 0.5*(sigmoid_D - sigmoid_D0) * d*d*sigmoid_s*sigmoid_s*sigmoid_s * std::exp (-d*sigmoid_s));

}
/// @brief Compute the derivative of depth- and bilayer-dependent dielectric
core::Real
ImplicitMembraneCoulomb::compute_deps_dr(
	core::Real const fi, /* fractional hydration of atom i */
	core::Real const fj, /* fractional hydration of atom j */
	core::Real const d
) const {

	// E = (1-f_ij)*d*D0 + f_ij*d*D
	// dE/dr = fD - fD0 + D0 = f(D-D0) + D0
	//E = (1-fij)*eps_mp + fij*eps_sol
	//dE/dr = (1-fij)eps_mp' + fij*eps_sol'

	core::Real f_ij = sqrt(fi*fj);
	core::Real D0_sol( 6.0 );
	core::Real D_sol( 80.0 );

	core::Real D0_mp( 3.0 );
	core::Real D_mp( 10.0 );
	core::Real slope( 0.40 );

	core::Real deps_dr = (1.0-f_ij)*sigmoidal_deps_dr(D0_mp, D_mp, slope, d);
	deps_dr = deps_dr + (f_ij)*sigmoidal_deps_dr(D0_sol, D_sol, slope, d);

	return deps_dr;
}

/// @brief Evaluate atom pair electrostatics energy (pass through)
core::Real
ImplicitMembraneCoulomb::eval_atom_atom_fa_elecE(
	core::Vector const & i_xyz,
	core::Real const i_charge,
	core::Real const i_hyd,
	core::Vector const & j_xyz,
	core::Real const j_charge,
	core::Real const j_hyd
) const {

	core::Real d2;
	return eval_atom_atom_fa_elecE(i_xyz, i_charge, i_hyd, j_xyz, j_charge, j_hyd, d2);
}

/// @brief Evaluate atom pair electrostatics energy
core::Real
ImplicitMembraneCoulomb::eval_atom_atom_fa_elecE(
	core::Vector const & i_xyz,
	core::Real const i_charge,
	core::Real const i_hyd,
	core::Vector const & j_xyz,
	core::Real const j_charge,
	core::Real const j_hyd,
	DistanceSquared & d2
) const {


	d2 = i_xyz.distance_squared( j_xyz );
	//have to substract the values when it is in water totally .
	if ( d2 > max_dis2_  ) {

		return 0.0;

	} else if ( d2 < low_poly_start2_ ) {

		return i_charge * j_charge * (compute_min_dis_score(i_hyd, j_hyd) - compute_min_dis_score(1.0,1.0));

	} else if ( d2 < low_poly_end2_ ) {

		core::Real score_sol(numeric::eval_cubic_polynomial( std::sqrt(d2), compute_lowpoly(1.0, 1.0) ));
		core::Real score_mp(numeric::eval_cubic_polynomial( std::sqrt(d2), compute_lowpoly(i_hyd, j_hyd) ));

		return i_charge * j_charge * (score_mp - score_sol);

	} else if ( d2 > hi_poly_start2_ ) {
		core::Real score_sol(numeric::eval_cubic_polynomial( std::sqrt(d2), compute_hipoly(1.0, 1.0) ));
		core::Real score_mp(numeric::eval_cubic_polynomial( std::sqrt(d2), compute_hipoly(i_hyd, j_hyd)));
		return i_charge * j_charge * (score_mp - score_sol);

	} else {
		core::Real C2(0.0), C2_sol(0.0);
		C2 = C1_/(max_dis_ * compute_depth_and_bilayer_dep_dielectric( i_hyd, j_hyd, max_dis_ ));
		C2_sol = C1_/(max_dis_ * compute_depth_and_bilayer_dep_dielectric( 1.0, 1.0, max_dis_ ));

		core::Real d = std::sqrt(d2);
		core::Real score_mp((C0_ /( d*compute_depth_and_bilayer_dep_dielectric( i_hyd, j_hyd, d ))) - C2);
		core::Real score_sol((C0_ /( d*compute_depth_and_bilayer_dep_dielectric( 1.0, 1.0, d ))) - C2_sol);
		return( i_charge*j_charge*(score_mp - score_sol) );
	}

}

/// @brief Get the key numeric value for derivative calculations
/// i.e. the derivative of energy with respect to distance divided by the distance
core::Real
ImplicitMembraneCoulomb::eval_dfa_elecE_dr_over_r(
	core::Real const dis2,
	core::Real const q1,
	core::Real const q2,
	core::Real const i_hyd,
	core::Real const j_hyd
) const {

	if ( dis2 > max_dis2_ ) return 0.0;
	else if ( dis2 < low_poly_start2_ ) return 0.0;

	core::Real q1q2 = q1*q2;

	if ( dis2 > low_poly_end2_ && dis2 < hi_poly_start2_ ) {
		core::Real d = std::sqrt(dis2);
		Real eps_d = compute_depth_and_bilayer_dep_dielectric( i_hyd, j_hyd, d );
		Real deps_d = compute_deps_dr( i_hyd, j_hyd, d );
		Real eps_sol = compute_depth_and_bilayer_dep_dielectric( 1.0, 1.0, d );
		Real deps_sol = compute_deps_dr( 1.0, 1.0, d );

		return dEfac_ * q1q2 *( ((eps_d + d*deps_d) / ( d*dis2*eps_d*eps_d )) - ((eps_sol + d*deps_sol) / ( d*dis2*eps_sol*eps_sol )) );

	} else if ( dis2 < low_poly_end2_ ) {
		Real d = std::sqrt( dis2 );
		return (numeric::cubic_polynomial_deriv( d, compute_lowpoly(i_hyd, j_hyd)) - numeric::cubic_polynomial_deriv( d, compute_lowpoly(1.0, 1.0))) * q1q2 / d;
	} else {
		Real d = std::sqrt( dis2 );
		return (numeric::cubic_polynomial_deriv( d, compute_hipoly(i_hyd, j_hyd) ) - numeric::cubic_polynomial_deriv( d, compute_hipoly(1.0, 1.0) ))* q1q2 / d;
	}
}

/// @brief Get the key numeric value for derivative calculations
/// i.e. the derivative of energy with respect to membrane depth
core::Vector
ImplicitMembraneCoulomb::eval_dfa_elecE_df(
	core::Real const dis2,
	core::Real const q1,
	core::Real const q2,
	core::Real const i_hyd,
	core::Real const j_hyd,
	core::Vector const di_hyd_df,
	core::Vector const dj_hyd_df
) const {

	if ( dis2 > max_dis2_ ) return core::Vector(0.0,0.0,0.0);
	else if ( dis2 < low_poly_start2_ ) return core::Vector(0.0,0.0,0.0);

	core::Real q1q2 = q1*q2;
	core::Real f_ij = sqrt( i_hyd*j_hyd );
	if ( dis2 > low_poly_end2_ && dis2 < hi_poly_start2_ ) {
		core::Real d = std::sqrt(dis2);
		core::Real eps_d = compute_depth_and_bilayer_dep_dielectric( i_hyd, j_hyd, d );
		core::Real eps_memb = compute_depth_and_bilayer_dep_dielectric( 0.0, 0.0, d );
		core::Real eps_sol = compute_depth_and_bilayer_dep_dielectric( 1.0, 1.0, d );


		return C0_ * q1q2 * ( ( eps_memb - eps_sol ) / ( d * eps_d * eps_d * f_ij ) ) * ( (i_hyd*dj_hyd_df) + (j_hyd*di_hyd_df) );

	} else if ( dis2 < low_poly_end2_ ) {
		Real d = std::sqrt( dis2 );
		core::Real score_mp(numeric::eval_cubic_polynomial( d, compute_lowpoly(i_hyd, j_hyd) ));
		return ( score_mp * q1q2/f_ij * ( (i_hyd*dj_hyd_df) + (j_hyd*di_hyd_df) ) );

	} else {
		Real d = std::sqrt( dis2 );
		core::Real score_mp(numeric::eval_cubic_polynomial( d, compute_hipoly(i_hyd, j_hyd) ));
		return ( score_mp * q1q2/f_ij * ( (i_hyd*dj_hyd_df) + (j_hyd*di_hyd_df) ) );
	}
}

} //energy_methods
} //core
