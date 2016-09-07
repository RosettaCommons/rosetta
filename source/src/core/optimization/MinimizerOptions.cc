// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/MinimizerOptions.cc
/// @brief  Minimizer options class implementation
/// @author Phil Bradley


// Unit headers
#include <core/optimization/MinimizerOptions.hh>

// Project headers
#include <basic/options/option.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {

/// @details Auto-generated virtual destructor
MinimizerOptions::~MinimizerOptions() = default;

MinimizerOptions::MinimizerOptions(
	std::string  min_type_in,
	Real const minimize_tolerance_in,
	bool const use_nblist_in,
	bool const deriv_check_in,
	bool const deriv_check_verbose_in
):
	max_iter_(2000),
	min_type_(std::move( min_type_in )),
	minimize_tolerance_( minimize_tolerance_in ),
	use_nblist_( use_nblist_in ),
	nblist_auto_update_( false ),
	deriv_check_( deriv_check_in ),
	deriv_check_verbose_( deriv_check_verbose_in ),
	deriv_check_to_stdout_( true ), // by default, send the result of a deriv-check call to the screen.
	// should make these configurable of course
	silent_( false ),
	gmax_cutoff_for_convergence_( 1.0 ),
	ax_init_( 0.0 ),
	xx_init_( 0.1 ),
	bx_init_( 0.2 ),
	brent_abs_tolerance_( 0.01 ),
	linmin_deriv_cutoff_( 0.0001 ),
	ga_mutation_probability_( 0.5 )
{
	using namespace basic::options;
	if ( option[ OptionKeys::run::nblist_autoupdate ].user() ) {
		nblist_auto_update_ = option[ OptionKeys::run::nblist_autoupdate ]();
	}
	if ( option[ OptionKeys::optimization::default_max_cycles ].user() ) {
		max_iter_ = option[ OptionKeys::optimization::default_max_cycles ]();
	}
}


MinimizerOptionsOP
MinimizerOptions::clone() const
{
	MinimizerOptionsOP minoptop( new MinimizerOptions(
		min_type_, minimize_tolerance_, use_nblist_,
		deriv_check_, deriv_check_verbose_ ) );
	if ( nblist_auto_update_ ) minoptop->nblist_auto_update_ = true;
	if (  deriv_check_to_stdout_ ) minoptop->deriv_check_to_stdout_ = true;

	return minoptop;
}

/////////////////////////////////////////////////////////////////////////////
// high-level params

// the min-type, eg "lbfgs_armijo_nonmonotone", "linmin"
std::string const &
MinimizerOptions::min_type() const
{
	return min_type_;
}

void
MinimizerOptions::min_type( std::string min_type_in )
{
	min_type_ = min_type_in;
}

std::string &
MinimizerOptions::min_type()
{
	return min_type_;
}

void
MinimizerOptions::deriv_check( bool deriv_check_in )
{
	deriv_check_ = deriv_check_in;
}

void
MinimizerOptions::deriv_check_to_stdout( bool setting )
{
	deriv_check_to_stdout_ = setting;
}


bool
MinimizerOptions::deriv_check() const
{
	return deriv_check_;
}


bool
MinimizerOptions::deriv_check_verbose() const
{
	return deriv_check_verbose_;
}

bool
MinimizerOptions::deriv_check_to_stdout() const
{
	return deriv_check_to_stdout_;
}


// the tolerance for dfpmin, dfpmin_atol etc
Real
MinimizerOptions::minimize_tolerance() const
{
	return minimize_tolerance_;
}

Real &
MinimizerOptions::minimize_tolerance()
{
	return minimize_tolerance_;
}

void
MinimizerOptions::minimize_tolerance( Real minimize_tolerance_in )
{
	minimize_tolerance_=minimize_tolerance_in;
}


bool
MinimizerOptions::use_nblist() const
{
	return use_nblist_;
}

void
MinimizerOptions::use_nblist( bool use_nblist_in )
{
	use_nblist_ = use_nblist_in;
}

bool
MinimizerOptions::nblist_auto_update() const {
	return nblist_auto_update_;
}

void
MinimizerOptions::nblist_auto_update( bool setting ) {
	nblist_auto_update_ = setting;
}

bool
MinimizerOptions::silent() const {
	return silent_;
}

void
MinimizerOptions::silent( bool setting ) {
	silent_ = setting;
}

/////////////////////////////////////////////////////////////////////////////
// low-level params

Real
MinimizerOptions::gmax_cutoff_for_convergence() const {
	return gmax_cutoff_for_convergence_;
}

void
MinimizerOptions::gmax_cutoff_for_convergence( Real setting ) {
	gmax_cutoff_for_convergence_ = setting;
}

// bracketing params for linmin
Real
MinimizerOptions::ax_init() const
{
	return ax_init_;
}

Real
MinimizerOptions::xx_init() const
{
	return xx_init_;
}

Real
MinimizerOptions::bx_init() const
{
	return bx_init_;
}

// abs tolerance for brent
Real
MinimizerOptions::brent_abs_tolerance() const
{
	return brent_abs_tolerance_;
}

/// @brief The derivative cutoff used for Brent.
///
Real
MinimizerOptions::linmin_deriv_cutoff() const {
	return linmin_deriv_cutoff_;
}

/// @brief Set the derivative cutoff used for Brent.
///
void
MinimizerOptions::linmin_deriv_cutoff(
	core::Real const &val
) {
	linmin_deriv_cutoff_ = val;
	return;
}


int MinimizerOptions::max_iter() const { return max_iter_; }
void MinimizerOptions::max_iter(int n) { max_iter_ = n; }


Real MinimizerOptions::ga_mutation_probability() const { return ga_mutation_probability_; }
void MinimizerOptions::ga_mutation_probability(Real p) { ga_mutation_probability_ = p; }


} // namespace optimization
} // namespace core
