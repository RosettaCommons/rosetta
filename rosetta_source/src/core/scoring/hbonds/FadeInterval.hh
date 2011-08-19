// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
/// class FadeInterval
/// stores an "fading interval" [a b c d] with a <= b <= c <= d
/// and computes the fraction of containment for x, which is defined to be
/// 0 if x is outside of (a,d), 1 if x is inside of [b,c],
/// and a linear ramp otherwise.
///    	___/-----\___
/// i.e. (x-a)/(b-a) for x in (a,b), and (d-x)/(d-c) for x in (c,d)
/// This is used to ensure that hbond scoring as a sum Er + ExH + ExD goes to zero at the edges.
///
/// Notes about discontinuities:
///   if x equals a, b, c, or d then deriv = 0
///   if x == a == b then value = 0
///   if x == c == d and a < c then value = 1
///   if x == c == d and a == c then value = 0
///   In particular if a == b == c == d then for all x, value == deriv == 0
////////////////////////////
class FadeInterval : public utility::pointer::ReferenceCount {

public:
	/// Constructor

	FadeInterval();

	FadeInterval(
		Real const min0,
		Real const fmin,
		Real const fmax,
		Real const max0);

	FadeInterval(
		std::string const & name,
		Real const min0,
		Real const fmin,
		Real const fmax,
		Real const max0);

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

};

} // hbonds
} // scoring
} // core

#endif
