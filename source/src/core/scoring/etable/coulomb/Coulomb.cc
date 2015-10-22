// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/coulomb/Coulomb.hh
/// @brief  Evaluate Coulombic potential
/// @author Phil Bradley, modifed by James Gleixner
/// @author Matthew O'Meara

// unit headers
#include <core/scoring/etable/coulomb/Coulomb.hh>

// project headers
#include <core/scoring/methods/EnergyMethodOptions.hh>

// numeric headers
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>

namespace core {
namespace scoring {
namespace etable {
namespace coulomb {


////////////////////////////////////////////////////////////////////////////
Coulomb::Coulomb( methods::EnergyMethodOptions const & options ):
	max_dis_( options.elec_max_dis() ),
	min_dis_( options.elec_min_dis() ),
	smooth_fa_elec_( options.smooth_fa_elec() ),
	die_( options.elec_die() ),
	no_dis_dep_die_( options.elec_no_dis_dep_die() ),
	sigmoidal_die_( options.elec_sigmoidal_die() )
{
	options.elec_sigmoidal_die_params(sigmoidal_D_, sigmoidal_D0_, sigmoidal_S_);
	initialize();
}


////////////////////////////////////////////////////////////////////////////
Coulomb::Coulomb( Coulomb const & src ): ReferenceCount(),
	max_dis_( src.max_dis_ ),
	min_dis_( src.min_dis_ ),
	smooth_fa_elec_( src.smooth_fa_elec_ ),
	die_( src.die_ ),
	no_dis_dep_die_( src.no_dis_dep_die_ ),
	sigmoidal_die_( src.sigmoidal_die_ ),
	sigmoidal_D_( src.sigmoidal_D_ ),
	sigmoidal_D0_( src.sigmoidal_D0_ ),
	sigmoidal_S_( src.sigmoidal_S_ )
{
	initialize();
}

CoulombOP
Coulomb::clone() const
{
	CoulombOP coulomb_ptr = CoulombOP( new Coulomb( *this ) );
	return coulomb_ptr;
}

void
Coulomb::initialize() {
	// Must have already initialized max_dis_, min_dis_, die_, and no_dis_dep_die_

	//max_dis_ = 5.5;
	max_dis2_ = max_dis_ * max_dis_;
	//min_dis_ = 1.5;
	min_dis2_ = min_dis_ * min_dis_ ;

	// default dielectric is 10r

	C0_ = 322.0637 ;
	if ( sigmoidal_die_ ) {
		C1_ = C0_ ;
		C2_ = C1_ / (max_dis_*sigmoid_eps (max_dis_));
		min_dis_score_ = C1_ / (min_dis_*sigmoid_eps (min_dis_)) - C2_ ;
		dEfac_ = -1.0 * C0_ ;
	} else if ( no_dis_dep_die_ ) {
		C1_ = C0_ / die_ ;
		C2_ = C1_ / max_dis_ ;
		min_dis_score_ = C1_ / min_dis_ - C2_ ;
		dEfac_ = -1.0 * C0_ / die_ ;
	} else {
		C1_ = C0_ / die_ ;
		C2_ = C1_ / max_dis2_ ;
		min_dis_score_ = C1_ / min_dis2_ - C2_ ;
		dEfac_ = -2.0 * C0_ / die_ ;
	}

	if ( smooth_fa_elec_ ) {
		low_poly_start_ = min_dis_ - 0.25;
		low_poly_end_   = min_dis_ + 0.25;
		low_poly_start2_ = low_poly_start_ * low_poly_start_;
		low_poly_end2_   = low_poly_end_ * low_poly_end_;

		// scope low polynomial
		{
			using namespace numeric::interpolation::spline;
			Real low_poly_end_score(0.0), low_poly_end_deriv(0.0);

			if ( sigmoidal_die_ ) {
				Real eps_low = sigmoid_eps (low_poly_end_);
				Real deps_low = sigmoid_deps_dr (low_poly_end_);
				low_poly_end_score = C1_ / (low_poly_end_*eps_low) - C2_;
				low_poly_end_deriv = -(C0_*(eps_low + low_poly_end_*deps_low))/(low_poly_end_*low_poly_end_*eps_low*eps_low);
			} else if ( no_dis_dep_die_ ) {
				low_poly_end_score = C1_ / low_poly_end_ - C2_;
				low_poly_end_deriv = -1 * C1_ / low_poly_end2_ ;
			} else {
				low_poly_end_score = C1_ / low_poly_end2_ - C2_;
				low_poly_end_deriv = -2 * C1_ / ( low_poly_end2_ * low_poly_end_ );
			}
			SplineGenerator gen_low_poly(
				low_poly_start_, min_dis_score_, 0,
				low_poly_end_, low_poly_end_score, low_poly_end_deriv );
			InterpolatorOP interp_low( gen_low_poly.get_interpolator() );
			SimpleInterpolatorOP sinterp_low = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_low );
			if ( ! sinterp_low ) {
				utility_exit_with_message( "Hack Elec created non-simple-interpolator in initialize()" );
			}
			SplineParameters low_sp;
			low_sp.ylo  = sinterp_low->y()[ 1 ];
			low_sp.yhi  = sinterp_low->y()[ 2 ];
			low_sp.y2lo = sinterp_low->ddy()[ 1 ];
			low_sp.y2hi = sinterp_low->ddy()[ 2 ];
			low_poly_ = Etable::cubic_polynomial_from_spline( low_poly_start_, low_poly_end_, low_sp );
		}

		hi_poly_start_    = max_dis_ - 1.0;
		hi_poly_end_      = max_dis_;
		hi_poly_start2_   = hi_poly_start_ * hi_poly_start_;
		hi_poly_end2_     = hi_poly_end_ * hi_poly_end_;

		// scope hi polynomial
		{
			using namespace numeric::interpolation::spline;
			Real hi_poly_start_score(0.0), hi_poly_start_deriv(0.0);
			if ( sigmoidal_die_ ) {
				Real eps_hi = sigmoid_eps (hi_poly_start_);
				Real deps_hi = sigmoid_deps_dr (hi_poly_start_);
				hi_poly_start_score = C1_ / (hi_poly_start_*eps_hi) - C2_;
				hi_poly_start_deriv = -(C0_*(eps_hi + hi_poly_start_*deps_hi))/(hi_poly_start_*hi_poly_start_*eps_hi*eps_hi);
			} else if ( no_dis_dep_die_ ) {
				hi_poly_start_score = C1_ / hi_poly_start_ - C2_;
				hi_poly_start_deriv = -1 * C1_ / hi_poly_start2_ ;

			} else {
				hi_poly_start_score = C1_ / hi_poly_start2_ - C2_;
				hi_poly_start_deriv = -2 * C1_ / ( hi_poly_start2_ * hi_poly_start_ );
			}

			SplineGenerator gen_hi_poly(
				hi_poly_start_, hi_poly_start_score, hi_poly_start_deriv,
				hi_poly_end_, 0, 0 );
			InterpolatorOP interp_hi( gen_hi_poly.get_interpolator() );
			SimpleInterpolatorOP sinterp_hi = utility::pointer::dynamic_pointer_cast< numeric::interpolation::spline::SimpleInterpolator > ( interp_hi );
			if ( ! sinterp_hi ) {
				utility_exit_with_message( "Hack Elec created non-simple-interpolator in initialize()" );
			}
			SplineParameters hi_sp;
			hi_sp.ylo  = sinterp_hi->y()[ 1 ];
			hi_sp.yhi  = sinterp_hi->y()[ 2 ];
			hi_sp.y2lo = sinterp_hi->ddy()[ 1 ];
			hi_sp.y2hi = sinterp_hi->ddy()[ 2 ];
			hi_poly_ = Etable::cubic_polynomial_from_spline( hi_poly_start_, hi_poly_end_, hi_sp );
		}
	} else {
		low_poly_start_ = min_dis_;     low_poly_start2_ = std::pow( low_poly_start_, 2 );
		low_poly_end_   = min_dis_ / 2; low_poly_end2_   = std::pow( low_poly_end_, 2 );
		hi_poly_start_  = max_dis_;     hi_poly_start2_  = std::pow( hi_poly_start_, 2 );
	}
}

}
}
}
}
