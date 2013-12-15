// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/MixtureFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


#include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/MixtureFunc.hh>

#include <core/types.hh>
#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

// C++ Headers

#include <iostream>
#include <algorithm>


// option key includes

#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace func {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::scoring::constraints;

	void
	MixtureFunc::read_data( std::istream & in ) {
		in 	>> anchor_ >> gaussian_param_ >> exp_param_ >> mixture_param_ >> bg_mean_ >> bg_sd_;
		verify_parameters_();
	}

	void MixtureFunc::verify_parameters_() {
		// quick check in case bg_mean_ and bg_sd_ aren't properly defined!
		if ( bg_mean_ <= 0.001 ) { // should never be true, atoms aren't this close!
			bg_mean_ =  21.0477;
			bg_sd_   =  6.60232;
		}


		rmax_ = anchor_ + 8;
		fmax_ = func_(rmax_);
		// std::cout << "(" << rmax_ << "," << fmax_ << ")" << std::endl;
	}

	/// @brief private
	Real MixtureFunc::func_( Real local_x ) const {
		// halved because an exponential is a one-sided distribution.
		Real exp_score   = 0.5 * dexponential( local_x, anchor_, exp_param_, mixture_param_ );
		Real gauss_score = dgaussian   ( local_x, anchor_, gaussian_param_, (1.0 - mixture_param_) );
		Real bg_score    = dgaussian   ( local_x, bg_mean_, bg_sd_, 1.0 );

		Real score;
		// std::cout << "score(" << local_x << ") = -1 * log( (" << exp_score << " + " << gauss_score
		// 					<< ") / " << exp(bg_score) << ")" << std::endl;
		if ( option[ basic::options::OptionKeys::constraints::normalize_mixture_func ]() ) {
			// brutal hack, make potential approach zero as we go from 3.8 -> 3.5 angstroms ...
			static const Real upper( 3.8 );
			static const Real lower( 3.5 );
			if ( local_x < upper ) {
				Real factor = std::min( (upper - local_x) / (upper - lower), 1.0 );
				bg_score = factor * (exp_score + gauss_score);
			}
			score = -1 * log( (exp_score + gauss_score) / bg_score );
		} else {
			score = -1 * log  (exp_score + gauss_score);
			// another brutal hack: if log is negative, set it to zero. This
			// happens when the estimated probabilities are greater than 1.0,
			// which is a numerical artifact of the function fitting process.
			if ( score < 0 ) score = 0;
		}
		// std::cout << "score = " << score << std::endl << std::endl;
		return score;
	}

	Real
	MixtureFunc::func( Real const x ) const	{
		Real local_x = x;
		// Real exp_score   = dexponential( local_x, anchor_, exp_param_, mixture_param_ );
		// Real gauss_score = dgaussian   ( local_x, anchor_, gaussian_param_, (1 - mixture_param_) );
		// Real bg_score    = dgaussian( local_x, bg_mean_, bg_sd_, 1 );

		Real score = func_(local_x);
		if ( local_x < rmax_ ) {
			score = func_(local_x);
		} else {
			score = fmax_ + local_x - rmax_;
		}

		//std::cout << "score(" << local_x << ") = " << func_(local_x) << std::endl;;

		if ( option[ basic::options::OptionKeys::constraints::penalize_mixture_func ]() ) {
			return score;
		} else {
			return std::min( score, 0.0 );
		}
	} // MixtureFunc::func( Real const x )


	Real
	MixtureFunc::dfunc( Real const x ) const {
		//Real g, h, g_prime, h_prime;
		//Real df = dfunc_component( x, g, h, g_prime, h_prime );

		//Real local_x = x;
		Real df = estimate_dfunc( x );
		// std::cout << "df = " << df << std::endl;

		return df;
	} // dfunc_component

	Real
	MixtureFunc::dfunc_component(
		Real const x,
		Real & g,
		Real & h,
		Real & g_prime,
		Real & h_prime
	)	const	{
		Real exp_deriv   = exponential_deriv( x, anchor_, exp_param_, mixture_param_ );
		Real gauss_deriv = gaussian_deriv   ( x, anchor_, gaussian_param_, (1 - mixture_param_) );

		Real exp_score   = dexponential( x, anchor_, exp_param_, mixture_param_ );
		Real gauss_score = dgaussian   ( x, anchor_, gaussian_param_, (1 - mixture_param_) );

		// f(x)  = log(g(x)) - log(h(x))
		// f'(x) = g'(x) / g(x) - h'(x) / h(x)
		g = exp_score + gauss_score;
		h = dgaussian( x, bg_mean_, bg_sd_, 1 );
		g_prime = exp_deriv + gauss_deriv;
		h_prime = gaussian_deriv( x, bg_mean_, bg_sd_, 1 );
		Real dfunc = g_prime / g - h_prime / h;

		return dfunc;
	} // dfunc_component

	void MixtureFunc::show( std::ostream& out ) const {
		using namespace ObjexxFCL::format;

		float start = 2;
		float end   = 20;
		float res   = 0.1;
		int width   = 16;
		int precise = 4;
		out << A( width, "r" )
				<< A( width, "func"  )
				<< A( width, "dfunc" )
				<< A( width, "g" )
				<< A( width, "h" )
				<< A( width, "g_prime" )
				<< A( width, "h_prime" )
				<< A( width, "dfunc_est" )
				<< std::endl;
		for ( Real r = start; r <= end; r += res ) {
			Real g, h, g_prime, h_prime;
			Real dfunc = dfunc_component( r, g, h, g_prime, h_prime );

			out << I( width, r )
					<< F( width, precise, func(r)  )
					<< F( width, precise, dfunc )
					<< F( width, precise, g )
					<< F( width, precise, h )
					<< F( width, precise, g_prime )
					<< F( width, precise, h_prime )
					<< F( width, precise, estimate_dfunc(r) )
					<< std::endl;
		} // for ( Real r = start; r <= end; r += res )
	} // void show( std::ostream& out )

	void MixtureFunc::show_definition( std::ostream &out ) const {
		out << "MIXTUREFUNC " << anchor_ << ' ' << gaussian_param_
				<< ' ' << exp_param_ << ' ' << mixture_param_ << ' ' << bg_mean_
				<< ' ' << bg_sd_ << std::endl;
	} // show_definition

	Real MixtureFunc::calc_kl_divergence() const {
		float res   = 0.1;
		float start = 2;
		float end   = 20;

		Real kl_divergence = 0;
		for ( float r = start; r <= end; r += res ) {
			Real g, h, g_prime, h_prime;
			dfunc_component( r, g, h, g_prime, h_prime );
			kl_divergence += g * ( std::log(g) - std::log(h) );
		} // for ( float r = start; r <= end; r += res )

		return kl_divergence;
	} // calc_kl_divergence

} // namespace constraints
} // namespace scoring
} // namespace core
