// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/FadeInterval.cc
/// @brief  FadeInterval class, used a simplified cross term for the hbond score term
/// @author Jack Snoeyink
/// @author Matthew O'Meara

// Unit Headers
#include <core/scoring/hbonds/FadeInterval.hh>

// Project Headers
#include <core/types.hh>

// C++ Headers
#include <ostream>
#include <string>

#include <utility/assert.hh>

namespace core {
namespace scoring {
namespace hbonds {

/// @details Auto-generated virtual destructor
FadeInterval::~FadeInterval() {}

using std::string;
using std::ostream;

// For the default construction always return weight=0, deriv=0
// Proof: In the functions 'value_deriv' and 'value'
// if x < 0 then x < fmax_ and x <= min0_ so value and deriv both equal 0
// otherwise if x >= 0 then x >= max0_ so value and deriv both equal 0
FadeInterval::FadeInterval():
	name_("fade_zero"),
	min0_(0.0),
	fmin_(0.0),
	fmax_(0.0),
	max0_(0.0),
	dfade_min_(0.0),
	dfade_max_(0.0),
	smooth_(false)
{}

FadeInterval::FadeInterval(
	Real const min0,
	Real const fmin,
	Real const fmax,
	Real const max0,
	bool const smooth) :
	name_("UNKNOWN"),
	min0_(min0),
	fmin_(fmin),
	fmax_(fmax),
	max0_(max0),
	dfade_min_(1.0/static_cast<double>(fmin-min0)),
	dfade_max_(1.0/static_cast<double>(max0-fmax)),
	smooth_(smooth)
{
	debug_assert(min0 <= fmin && fmin <= fmax && fmax <= max0);
}

FadeInterval::FadeInterval(
	string const & name,
	Real const min0,
	Real const fmin,
	Real const fmax,
	Real const max0,
	bool const smooth) :
	name_(name),
	min0_(min0),
	fmin_(fmin),
	fmax_(fmax),
	max0_(max0),
	dfade_min_(1.0/static_cast<double>(fmin-min0)),
	dfade_max_(1.0/static_cast<double>(max0-fmax)),
	smooth_(smooth)
{
	debug_assert(min0 <= fmin && fmin <= fmax && fmax <= max0);
}

void
FadeInterval::value_deriv(
	Real const x,
	double &val,
	double &deriv) const
{
	//JSS  5 intervals --a-b---c-d--
	if ( x <= fmax_ ) {
		if ( x <= min0_ ) { val = deriv = 0.0; return; }       // in (-\infty, min0]
		if ( x >= fmin_ ) { val = 1.0; deriv = 0.0; return; }  // in [fmin, fmax]
		if ( smooth_ ) {
			double const z((x - min0_)* dfade_min_);
			val = z*z*(3-2*z);
			deriv = -6*z*(z-1)*dfade_min_;
		} else {
			deriv = dfade_min_; val = (x - min0_) * dfade_min_;  // in (min0,fmin)
		}
	} else {
		if ( x >= max0_ ) { val = deriv = 0.0; return; }       // in [max0, \infty)
		if ( smooth_ ) {
			double const z((x - fmax_) * dfade_max_);
			val = z*z*(2*z-3) + 1;
			deriv = 6*z*(z-1)*dfade_max_;
		} else {
			deriv = -dfade_max_; val = (max0_ - x) * dfade_max_; // in (fmax,max0)
		}
	}
}

double
FadeInterval::value(
	Real const x) const
{
	//JSS  5 intervals --a-b---c-d--
	if ( x <= fmax_ ) {
		if ( x <= min0_ ) return 0.0;      // in (-\infty, min0]
		if ( x >= fmin_ ) return 1.0;      // in [fmin, fmax]
		if ( smooth_ ) {
			double const z((x - min0_)* dfade_min_);
			return z*z*(3.0-2.0*z);
		} else {
			return (x - min0_) * dfade_min_; // in (min0,fmin)
		}
	} else {
		if ( x >= max0_ ) return 0.0;      // in [max0, \infty)
		if ( smooth_ ) {
			double const z((x - fmax_) * dfade_max_);
			return z*z*(2.0*z-3.0) + 1.0;
		} else {
			return (max0_ - x) * dfade_max_; // in (fmax,max0)
		}
	}
}

string
FadeInterval::get_name() const
{
	return name_;
}

Real
FadeInterval::get_min0() const
{
	return min0_;
}

Real
FadeInterval::get_fmin() const
{
	return fmin_;
}

Real
FadeInterval::get_fmax() const
{
	return fmax_;
}

Real
FadeInterval::get_max0() const
{
	return max0_;
}

bool
FadeInterval::get_smooth() const
{
	return smooth_;
}

bool
operator==(
	FadeInterval const & a,
	FadeInterval const & b)
{
	return
		a.min0_ == b.min0_ &&
		a.fmin_ == b.fmin_ &&
		a.fmax_ == b.fmax_ &&
		a.max0_ == b.max0_ &&
		a.smooth_ == b.smooth_;
}

bool
operator!=(
	FadeInterval const & a,
	FadeInterval const & b)
{
	return !(a == b);
}

ostream &
operator<< (
	ostream & out,
	FadeInterval const & fade_interval)
{
	fade_interval.show( out );
	return out;
}

void
FadeInterval::show(
	ostream & out ) const
{
	out << "(" << name_ << "; " << min0_ << ", " << fmin_ << ", " << fmax_ << ", " << max0_ << ")" << (smooth_ ? " smoothed" : "");
}

} // hbonds
} // scoring
} // core
