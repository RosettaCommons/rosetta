// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file core/scoring/hbonds/polynomial.cc
/// @brief Classes for hydrogen bond polynomial evaluation functions
/// @author Matthew O'Meara (mattjomeara@gmail.com


// Unit Headers
#include <core/scoring/hbonds/polynomial.hh>

// Package headers
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/hbonds/types.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/conversions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>
#include <string>
#include <cmath>

namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt;

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
	polynomial_name_(polynomial_name),
	geometric_dimension_(geometric_dimension),
	xmin_(xmin), xmax_(xmax), min_val_(min_val), max_val_(max_val), root1_(root1), root2_(root2),
	degree_(degree),
	coefficients_(coefficients)
{}

Polynomial_1d::Polynomial_1d(Polynomial_1d const & src):
	utility::pointer::ReferenceCount( src ),
	polynomial_name_(src.polynomial_name_),
	geometric_dimension_(src.geometric_dimension_),
	xmin_(src.xmin_), xmax_(src.xmax_), root1_(src.root1_), root2_(src.root2_),
	degree_(src.degree_),
	coefficients_(src.coefficients_)
{}

Polynomial_1d::~Polynomial_1d(){}

HBGeoDimType
Polynomial_1d::geometric_dimension() const
{
	return geometric_dimension_;
}

string
Polynomial_1d::name() const
{
	return polynomial_name_;
}

Real
Polynomial_1d::xmin() const
{
	return xmin_;
}

Real
Polynomial_1d::xmax() const
{
	return xmax_;
}

Real
Polynomial_1d::min_val() const
{
	return min_val_;
}

Real
Polynomial_1d::max_val() const
{
	return max_val_;
}

Real
Polynomial_1d::root1() const
{
	return root1_;
}

Real
Polynomial_1d::root2() const
{
	return root2_;
}

Size
Polynomial_1d::degree() const
{
	return degree_;
}

vector1< Real > const &
Polynomial_1d::coefficients() const
{
	return coefficients_;
}

////////////////////////////////////////////////////////////////////////////////
/// @begin operator()
///
/// @brief evaluate the polynomial and its derivative.
///
/// @detailed
///
/// @param  variable - [in] - evaluate polynomial(value)
/// @param  value - [out] - returned output
/// @param  deriv - [out] - returned output
///
/// @global_read
///
/// @global_write
///
/// @remarks
///  Note the coefficients must be in reverse order: low to high
///
///  Polynomial value and derivative using Horner's rule
///  value = Sum_(i = 1,...,N) [ coeff_i * variable^(i-1) ]
///  deriv = Sum_(i = 2,...,N) [ ( i - 1 ) * coeff_i * variable^(i-2) ]
///  JSS: Horner's rule for evaluating polynomials is based on rewriting the polynomial as:
///  JSS: p(x)  = a0 + x*(a1 + x*(a2 + x*(...  x*(aN)...)))
///  JSS: or value_k = a_k + x*value_k+1 for k = N-1 to 0
///  JSS: and the derivative is
///  JSS: deriv_k = value_k+1 + deriv_k+1 for k = N-1 to 1
///
/// @references
///
/// @authors Jack Snoeyink
/// @authors Matthew O'Meara
///
/// @last_modified Matthew O'Meara
/////////////////////////////////////////////////////////////////////////////////
void
Polynomial_1d::operator()(
	double const variable,
	double & value,
	double & deriv) const
{
	if(variable <= xmin_){
		value = min_val_;
		deriv = 0.0;
		return;
	}
	if(variable >= xmax_){
		value = max_val_;
		deriv = 0.0;
		return;
	}

	value = coefficients_[1];
	deriv = 0.0;
	for(Size i=2; i <= degree_; i++){
		(deriv *= variable) += value;
		(value *= variable) += coefficients_[i];
	}
}

ostream &
operator<< ( ostream & out, const Polynomial_1d & poly ){
	poly.show( out );
	return out;
}

void
Polynomial_1d::show( ostream & out ) const{
	out << polynomial_name_ << " "
			<< "geometric dimension:" <<HBondTypeManager::name_from_geo_dim_type(geometric_dimension_) << " "
			<< "domain:(" << xmin_ << "," << xmax_ << ") "
			<< "out_of_range_vals:(" << min_val_ << "," << max_val_ << ") "
			<< "roots:[" << root1_ << "," << root2_ << "] "
			<< "degree:" << degree_ << " "
			<< "y=";
	for(Size i=1; i <= degree_; ++i){
		if (i >1){
			if (coefficients_[i] > 0 ){
				out << "+";
			} else if (coefficients_[i] < 0 ){
				out << "-";
			} else{
				continue;
			}
		}
		out << std::abs(coefficients_[i]);
		if (degree_-i >1){
			out << "x^" << degree_-i;
		} else if (degree_-i == 1){
			out << "x";
		} else {}
	}
}

std::string
Polynomial_1d::show_values( ) const{

	using numeric::conversions::radians;

	Real energy, deriv;
	stringstream out;
	out << polynomial_name_;

	switch(geometric_dimension_){
	case(hbgd_NONE):
		out << " NONE ";
		break;
	case(hbgd_AHdist):
		out << " AHdist ";
		for ( Size i=120; i<= 300; i+=2 ){
			Real const AHdist( 0.01 * i );
			operator()(AHdist, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	case(hbgd_cosBAH):
		out << " cosBAH ";
		for ( Size i=60; i<= 180; i+=2 ) {
			Real const xH( cos( radians( 180.0 - i ) ) );
			operator()(xH, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	case(hbgd_cosAHD):
		for ( Size i=60; i<= 180; i+=2 ) {
			Real const xD( cos( radians( 180.0 - i ) ) );
			operator()(xD, energy, deriv);
			out << F(9, 3, energy) << " ";
		}
		break;
	case(hbgd_chi):
		// just a guess, needs checking !
		for ( Size i=0; i<=360; i+=6) {
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
