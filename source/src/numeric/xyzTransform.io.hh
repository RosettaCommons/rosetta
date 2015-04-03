// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzTransform.io.hh
/// @brief  xyzTransform input/output functions
/// @author will sheffler


#ifndef INCLUDED_numeric_xyzTransform_io_hh
#define INCLUDED_numeric_xyzTransform_io_hh


// Package headers
#include <numeric/xyzTransform.hh>
#include <numeric/xyzMatrix.io.hh>
#include <numeric/xyzVector.io.hh>
#include <ObjexxFCL/format.hh>
// C++ headers
#include <iostream>


namespace numeric {


/// @brief stream << xyzTransform output operator
template< typename T >
std::ostream &
operator <<( std::ostream & stream, xyzTransform< T > const & m )
{
	using namespace ObjexxFCL::format;
	T angle;
	xyzVector<T> axis = rotation_axis(m.R,angle);
	// if(axis.x()<0) { axis = -axis; angle = numeric::constants::d::pi_2-angle; }
	xyzVector<T> dirn = m.t.length() > 0.000001 ? m.t.normalized() : xyzVector<T>(0,0,0);
	stream << "axis: ( " << F(10,7,axis.x()) << " " << F(10,7,axis.y()) << " " << F(10,7,axis.z()) << " ) " << F(10,6,numeric::conversions::degrees(angle)) << " ";
	stream << "dirn: ( " << F(10,5,dirn.x()) << " " << F(10,5,dirn.y()) << " " << F(10,5,dirn.z()) << " ) * " << F(10,6,m.t.length());
	return stream;
}


/// @brief stream >> xyzTransform input operator
/// @note Reads row-ordered matrix elements from one or multiple lines
template< typename T >
std::istream &
operator >>( std::istream & stream, xyzTransform< T > & m )
{
	stream >> m.R >> m.t;
	return stream;
}


} // namespace numeric


#endif // INCLUDED_numeric_xyzTransform_io_HH
