// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/polynomial.hh
/// @brief  Polynomial evaluation class
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_numeric_polynomial_hh
#define INCLUDED_numeric_polynomial_hh

// Unit Headers
#include <numeric/polynomial.fwd.hh>

// Numeric Headers
#include <numeric/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <iostream>

namespace numeric {

class Polynomial_1d : public utility::pointer::ReferenceCount {

public:
	//Polynomial_1d();

	Polynomial_1d(
		std::string const & polynomial_name,
		Real const xmin,
		Real const xmax,
		Real const min_val,
		Real const max_val,
		Real const root1,
		Real const root2,
		Size degree,
		utility::vector1< Real > const & coefficients);

	Polynomial_1d(Polynomial_1d const & src );

	virtual ~Polynomial_1d();

	// check that the polynomial parameters are well formed
	void
	check_invariants() const;

	std::string
	name() const;

	Real
	xmin() const;

	Real
	xmax() const;

	Real
	min_val() const;

	Real
	max_val() const;

	Real
	root1() const;

	Real
	root2() const;

	Size
	degree() const;

	utility::vector1< Real > const &
	coefficients() const;

	///@brief Evaluate the polynomial and its derivative
	void
	operator()(
		double const variable,
		double & value,
		double & deriv) const;

	/// @brief just evaluate the polynomial, w/o derivatives
	double
	eval( double const variable );

	void
	show( std::ostream & out ) const;

private:
	std::string polynomial_name_;
	Real xmin_;
	Real xmax_;
	Real min_val_;
	Real max_val_;
	Real root1_;
	Real root2_;
	Size degree_;
	utility::vector1< Real > coefficients_;

};

std::ostream &
operator<< ( std::ostream & out, const Polynomial_1d & poly );

} // numeric

#endif // INCLUDED_numeric_polynomial_HH
