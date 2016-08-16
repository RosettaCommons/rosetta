// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/FadeInterval.hh
/// @brief  FadeInterval class, used a simplified cross term for the hbond score term
/// @author Jack Snoeyink
/// @author Matthew O'Meara

#ifndef INCLUDED_core_scoring_hbonds_FadeInterval_hh
#define INCLUDED_core_scoring_hbonds_FadeInterval_hh

// Unit headers
#include <core/scoring/hbonds/FadeInterval.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace core {
namespace scoring {
namespace hbonds {

/////////////////////////////
/// Classic FadeInterval
/////////////////////////////
/// stores an "fading interval" [a b c d] with a <= b <= c <= d
/// and computes the fraction of containment for x, which is defined to be
/// 0 if x is outside of (a,d), 1 if x is inside of [b,c],
/// and a linear ramp otherwise.
///     ___/-----\___
/// i.e. (x-a)/(b-a) for x in (a,b), and (d-x)/(d-c) for x in (c,d)
/// This is used to ensure that hbond scoring as a sum Er + ExH + ExD goes to zero at the edges.
///
/// Notes about discontinuities:
///   if x equals a, b, c, or d then deriv = 0
///   if x == a == b then value = 0
///   if x == c == d and a < c then value = 1
///   if x == c == d and a == c then value = 0
///   In particular if a == b == c == d then for all x, value == deriv == 0
///
//////////////////////////////
/// Smooth FadeInteval
//////////////////////////////
/// Rather than using a piecewise linear fading function,
/// use a piecewise sigmoid function to have a continuous derivative.
///
/// Look for a canonical sigmoid function f(x) such that,
///   f(0)  = 1    f(1)  = 0  // goes through the the knots
///   f'(0) = 0    f'(1) = 0  // is horizontal at the knots
///   a continuous differative
///   f(x-.5)-.5 is odd       // symmetric
///
/// I claim, f(x) = 2x^3 - 3x^2 + 1, satisfies these constraints:
///   f(0) = 2(0)^3 - 3(x)^2 + 1 = 1
///   f(1) = 2(1)^3 - 3(1)^2 + 1 = 2 - 3 + 1 = 0
///
///   f'(x) = 6x^2 - 6x = 6x(x-1)
///   f'(0) = 6(0)(0-1) = 0
///   f'(1) = 6(1)(1-1) = 0
///
///   a function g(x) is odd if g(-x) = -g(x)
///   f((-x)-.5)-.5 = 2(-x-.5)^3 - 3(-x-.5)^2 + 1 - .5
///                 = -2x^3 - 6x^2 - 4.5x - .5
///   -(f(x-.5)-.5) = -2(x-.5)^3 + 3(x-.5)^2 - 1 + .5
///                 = -2x^3 - 6x^2 - 4.5x - .5
///
/// Given the knots --a-b---c-d-- transform f(x) to fill a-b and c-d
/// to connect the linear regions.
///
///
///   a-b region:
///     let z(x) = (x-a)/(b-a) and z'(x) = 1/(b-a)
///     use g(x) = 1-f(-z(x)) = -2z^3 + 3z^2
///     and g'(x) = -6z(z-1)*z'(x)
///   c-d region:
///     let z(x) = (x-c)/(d-c) and z'(x) = 1/(d-c)
///     use g(x) = f(z) = 2z^3 - 3z^2 + 1
///     and g'(x) = 6z(z-1)*z'(x)
///
////////////////////////////
class FadeInterval : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FadeInterval();
	/// Constructor

	FadeInterval();

	FadeInterval(
		Real const min0,
		Real const fmin,
		Real const fmax,
		Real const max0,
		bool const smooth = false);

	FadeInterval(
		std::string const & name,
		Real const min0,
		Real const fmin,
		Real const fmax,
		Real const max0,
		bool const smooth = false);

	void
	value_deriv(
		Real const x,
		double &val,
		double &deriv
	) const;

	double
	value(
		Real const x
	) const;

	std::string
	get_name() const;

	Real
	get_min0() const;

	Real
	get_fmin() const;

	Real
	get_fmax() const;

	Real
	get_max0() const;

	bool
	get_smooth() const;

	friend
	bool
	operator==(FadeInterval const & a, FadeInterval const & b);

	friend
	bool
	operator!=(FadeInterval const & a, FadeInterval const & b);

	friend
	std::ostream &
	operator<< (std::ostream & out, FadeInterval const & fade_interval);

	void
	show( std::ostream & out ) const;

private:
	std::string const name_;
	Real const min0_;
	Real const fmin_;
	Real const fmax_;
	Real const max0_;
	double const dfade_min_;
	double const dfade_max_;
	bool const smooth_;

};

} // hbonds
} // scoring
} // core

#endif
