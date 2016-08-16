// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/polynomial.cc
/// @brief Classes for hydrogen bond polynomial evaluation functions
/// @author Matthew O'Meara (mattjomeara@gmail.com


// Unit Headers
#include <core/scoring/hbonds/polynomial.hh>

// Package headers
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/hbonds/types.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/polynomial.hh>
#include <numeric/conversions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>
#include <string>
#include <cmath>

namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format;

namespace core {
namespace scoring {
namespace hbonds {

using namespace std;
using namespace utility;

/// @brief ctor
Polynomial_1d::Polynomial_1d(
	string const & polynomial_name,
	HBGeoDimType const geometric_dimension,
	Real const xmin,
	Real const xmax,
	Real const min_val,
	Real const max_val,
	Real const root1,
	Real const root2,
	Size degree,
	vector1< Real > const & coefficients):
	numeric::Polynomial_1d(polynomial_name, xmin, xmax, min_val, max_val, root1, root2, degree, coefficients),
	geometric_dimension_(geometric_dimension)
{}

Polynomial_1d::Polynomial_1d(Polynomial_1d const & src):
	numeric::Polynomial_1d(src),
	geometric_dimension_(src.geometric_dimension_)
{}

HBGeoDimType
Polynomial_1d::geometric_dimension() const
{
	return geometric_dimension_;
}

ostream &
operator<< ( ostream & out, const Polynomial_1d & poly ){
	poly.show( out );
	return out;
}

void
Polynomial_1d::show( ostream & out ) const{
	out << name() << " "
		<< "geometric dimension:" << HBondTypeManager::name_from_geo_dim_type(geometric_dimension()) << " "
		<< "domain:(" << xmin() << "," << xmax() << ") "
		<< "out_of_range_vals:(" << min_val() << "," << max_val() << ") "
		<< "roots:[" << root1() << "," << root2() << "] "
		<< "degree:" << degree() << " "
		<< "y=";
	for ( Size i=1; i <= degree(); ++i ) {
		if ( i >1 ) {
			if ( coefficients()[i] > 0 ) {
				out << "+";
			} else if ( coefficients()[i] < 0 ) {
				out << "-";
			} else {
				continue;
			}
		}
		out << std::abs(coefficients()[i]);
		if ( degree()-i >1 ) {
			out << "x^" << degree()-i;
		} else if ( degree()-i == 1 ) {
			out << "x";
		} else { }
	}
}

std::string
Polynomial_1d::show_values( ) const{

	using numeric::conversions::radians;

	Real energy, deriv;
	stringstream out;
	out << name();

	switch(geometric_dimension_){
	case(hbgd_NONE) :
		out << " NONE ";
		break;
	case(hbgd_AHdist) :
		out << " AHdist ";
		for ( Size i=120; i<= 300; i+=2 ) {
			Real const AHdist( 0.01 * i );
			operator()(AHdist, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	case(hbgd_cosBAH) :
		out << " cosBAH ";
		for ( Size i=60; i<= 180; i+=2 ) {
			Real const xH( cos( radians( 180.0 - i ) ) );
			operator()(xH, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	case(hbgd_cosAHD) :
		//out << " cosAHD ";
		for ( Size i=60; i<= 180; i+=2 ) {
			Real const xD( cos( radians( 180.0 - i ) ) );
			operator()(xD, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	case(hbgd_AHD) :
		//out << " AHD ";
		// I am not sure how this should be implemented; I added this case block to remove compiler warning:
		// "enum value hbgd_AHD not handled in switch statement"
		// ~Labonte
		break;
	case(hbgd_chi) :
		// just a guess, needs checking !
		for ( Size i=0; i<=360; i+=6 ) {
			Real const chi( cos( radians( 180.0 - i ) ) );
			operator()(chi, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	}
	return out.str();
}

} // hbonds
} // scoring
} // core
