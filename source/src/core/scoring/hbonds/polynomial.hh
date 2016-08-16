// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/polynomial.hh
/// @brief  Polynomial objects for hydrogen bond score term
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_scoring_hbonds_polynomial_hh
#define INCLUDED_core_scoring_hbonds_polynomial_hh

// Unit Headers
#include <core/scoring/hbonds/polynomial.fwd.hh>
#include <numeric/polynomial.hh>
#include <core/scoring/hbonds/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <iostream>

namespace core {
namespace scoring {
namespace hbonds {

class Polynomial_1d : public numeric::Polynomial_1d {

public:
	//Polynomial_1d();

	Polynomial_1d(
		std::string const & polynomial_name,
		HBGeoDimType const geometric_dimension,
		Real const xmin,
		Real const xmax,
		Real const min_val,
		Real const max_val,
		Real const root1,
		Real const root2,
		Size degree,
		utility::vector1< Real > const & coefficients);

	Polynomial_1d(Polynomial_1d const & src);

	HBGeoDimType
	geometric_dimension() const;

	void
	show( std::ostream & out ) const;

	std::string
	show_values() const;

private:
	HBGeoDimType geometric_dimension_;
};

std::ostream &
operator<< ( std::ostream & out, const Polynomial_1d & poly );


} // hbonds
} // scoring
} // core

#endif // INCLUDED_core_scoring_hbonds_polynomial_HH
