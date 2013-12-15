// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/CircularSigmoidalFunc.cc
/// @brief Definition for functions used in definition of constraints.
/// @author Robert Vernon

#include <core/scoring/func/CircularSigmoidalFunc.hh>
#include <numeric/angle.functions.hh>
#include <core/types.hh>
#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif
#include <cmath>

#include <iostream>

namespace core {
namespace scoring {
namespace func {

Real
CircularSigmoidalFunc::func( Real const x ) const {
	Real const x0 = numeric::nearest_angle_radians(x,xC_)-xC_;

	//std::cout << "SIG_FUNC " << M_E << " " << x << " " << xC_ << " " << x0 << " " << o1_ << " " << o2_ << " " << m_ << " " << offset_ << " " << (1/(1+ std::pow( M_E, (-m_*(x0-o1_)) ))) << " " << (1/(1+ std::pow( M_E, (-m_*(x0-o2_)) ))) << std::endl;

	Real const z = offset_ + (1/(1+ std::pow( M_E, (-m_*(x0-o1_)) ))) - (1/(1+ std::pow( M_E, (-m_*(x0-o2_)) )));

	return z;
}

Real
CircularSigmoidalFunc::dfunc( Real const x ) const {

	Real const x0 = numeric::nearest_angle_radians(x,xC_)-xC_;

	Real const z = m_*std::pow(M_E,m_*x0+m_*o1_)/(2*std::pow(M_E,m_*x0+m_*o1_)+std::pow(M_E,2*m_*x0)+std::pow(M_E,2*m_*o1_))
	           	 - m_*std::pow(M_E,m_*x0+m_*o2_)/(2*std::pow(M_E,m_*x0+m_*o2_)+std::pow(M_E,2*m_*x0)+std::pow(M_E,2*m_*o2_));

	return z;
}

void
CircularSigmoidalFunc::read_data( std::istream & in ) {
	in >> xC_ >> m_ >> o1_ >> o2_;
}

void CircularSigmoidalFunc::show_definition( std::ostream & out ) const
{
	out << "CircularSigmoidalFunc " << xC_ << ' ' << m_ << ' ' << o1_ << ' ' << o2_;
}

} // namespace constraints
} // namespace scoring
} // namespace core
