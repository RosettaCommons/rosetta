// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/polynomial.cc
/// @brief Classes for polynomial evaluation functions
/// @author Matthew O'Meara (mattjomeara@gmail.com


// Unit Headers
#include <numeric/polynomial.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/types.hh>

// C++ headers
#include <iostream>
#include <string>
#include <cmath>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace numeric {

using std::string;
using std::ostream;
using utility::vector1;

/// @brief ctor
Polynomial_1d::Polynomial_1d(
	string const & polynomial_name,
	Real const xmin,
	Real const xmax,
	Real const min_val,
	Real const max_val,
	Real const root1,
	Real const root2,
	Size degree,
	vector1< Real > const & coefficients):
	polynomial_name_(polynomial_name),
	xmin_(xmin), xmax_(xmax), min_val_(min_val), max_val_(max_val), root1_(root1), root2_(root2),
	degree_(degree),
	coefficients_(coefficients)
{
	check_invariants();
}

Polynomial_1d::Polynomial_1d(Polynomial_1d const & src):
	utility::pointer::ReferenceCount( src ),
	polynomial_name_(src.polynomial_name_),
	xmin_(src.xmin_), xmax_(src.xmax_), root1_(src.root1_), root2_(src.root2_),
	degree_(src.degree_),
	coefficients_(src.coefficients_)
{
	check_invariants();
}

Polynomial_1d::~Polynomial_1d(){}

void
Polynomial_1d::check_invariants() const
{
	if ( xmin_ > xmax_ ) {
		std::stringstream msg;
		msg << "Polnomial_1d is badly formed because (xmin: '" << xmin_ << "') > (xmax: '" << xmax_ << "')";
		throw utility::excn::EXCN_Msg_Exception(msg.str() );
	}

	if ( coefficients_.size() == 0 ) {
		std::stringstream msg;
		msg << "Polnomial_1d is badly formed because no coefficients were provided";
		throw utility::excn::EXCN_Msg_Exception(msg.str() );
	}

	if ( coefficients_.size() != degree_ ) {
		std::stringstream msg;
		msg << "Polnomial_1d is badly formed because the degree was given to be '" << degree_ << "' while '" << coefficients_.size() << "' coefficients were provided, while they should be equal.";
		throw utility::excn::EXCN_Msg_Exception(msg.str() );
	}


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
///
/// @brief evaluate the polynomial and its derivative.
///
/// @details
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
/// @author Jack Snoeyink
/// @author Matthew O'Meara
///
/////////////////////////////////////////////////////////////////////////////////
void
Polynomial_1d::operator()(
	double const variable,
	double & value,
	double & deriv) const
{
	if ( variable <= xmin_ ) {
		value = min_val_;
		deriv = 0.0;
		return;
	}
	if ( variable >= xmax_ ) {
		value = max_val_;
		deriv = 0.0;
		return;
	}

	value = coefficients_[1];
	deriv = 0.0;
	for ( Size i=2; i <= degree_; i++ ) {
		(deriv *= variable) += value;
		(value *= variable) += coefficients_[i];
	}
}

double
Polynomial_1d::eval( double const variable )
{
	if ( variable <= xmin_ ) {
		return min_val_;
	}
	if ( variable >= xmax_ ) {
		return max_val_;
	}
	Real value = coefficients_[1];
	for ( Size i=2; i <= degree_; ++i ) {
		( value *= variable ) += coefficients_[i];
	}
	return value;
}


ostream &
operator<< ( ostream & out, const Polynomial_1d & poly ){
	poly.show( out );
	return out;
}

void
Polynomial_1d::show( ostream & out ) const{
	out << polynomial_name_ << " "
		<< "domain:(" << xmin_ << "," << xmax_ << ") "
		<< "out_of_range_vals:(" << min_val_ << "," << max_val_ << ") "
		<< "roots:[" << root1_ << "," << root2_ << "] "
		<< "degree:" << degree_ << " "
		<< "y=";
	for ( Size i=1; i <= degree_; ++i ) {
		if ( i >1 ) {
			if ( coefficients_[i] > 0 ) {
				out << "+";
			} else if ( coefficients_[i] < 0 ) {
				out << "-";
			} else {
				continue;
			}
		}
		out << std::abs(coefficients_[i]);
		if ( degree_-i >1 ) {
			out << "x^" << degree_-i;
		} else if ( degree_-i == 1 ) {
			out << "x";
		} else { }
	}
}

} // namespace

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
numeric::Polynomial_1d::Polynomial_1d() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::Polynomial_1d::save( Archive & arc ) const {
	arc( CEREAL_NVP( polynomial_name_ ) ); // std::string
	arc( CEREAL_NVP( xmin_ ) ); // Real
	arc( CEREAL_NVP( xmax_ ) ); // Real
	arc( CEREAL_NVP( min_val_ ) ); // Real
	arc( CEREAL_NVP( max_val_ ) ); // Real
	arc( CEREAL_NVP( root1_ ) ); // Real
	arc( CEREAL_NVP( root2_ ) ); // Real
	arc( CEREAL_NVP( degree_ ) ); // Size
	arc( CEREAL_NVP( coefficients_ ) ); // utility::vector1<Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::Polynomial_1d::load( Archive & arc ) {
	arc( polynomial_name_ ); // std::string
	arc( xmin_ ); // Real
	arc( xmax_ ); // Real
	arc( min_val_ ); // Real
	arc( max_val_ ); // Real
	arc( root1_ ); // Real
	arc( root2_ ); // Real
	arc( degree_ ); // Size
	arc( coefficients_ ); // utility::vector1<Real>
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::Polynomial_1d );
CEREAL_REGISTER_TYPE( numeric::Polynomial_1d )

CEREAL_REGISTER_DYNAMIC_INIT( numeric_polynomial )
#endif // SERIALIZATION
