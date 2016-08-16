// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyz.functions.fwd.hh
/// @brief  numeric::xyzVector and numeric::xyzMatrix functions forward declarations
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_xyz_functions_fwd_hh
#define INCLUDED_numeric_xyz_functions_fwd_hh


// Package headers
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/sphericalVector.fwd.hh>


namespace numeric {


// Forward
template< typename T > xyzVector< T > operator *( xyzMatrix< T > const & m, xyzVector< T > const & v );
template< typename T > xyzVector< T > product( xyzMatrix< T > const & m, xyzVector< T > const & v );
template< typename T > xyzVector< T > & inplace_product( xyzMatrix< T > const & m, xyzVector< T > & v );
template< typename T > xyzVector< T > transpose_product( xyzMatrix< T > const & m, xyzVector< T > const & v );
template< typename T > xyzVector< T > & inplace_transpose_product( xyzMatrix< T > const & m, xyzVector< T > & v );
template< typename T > xyzMatrix< T > outer_product( xyzVector< T > const & a, xyzVector< T > const & b );
template< typename T > xyzMatrix< T > inverse( xyzMatrix< T > const & a );
template< typename T > xyzMatrix< T > projection_matrix( xyzVector< T > const & v );
template< typename T > void dihedral_radians( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & p3, xyzVector< T > const & p4, T & angle );
template< typename T > T dihedral_radians( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & p3, xyzVector< T > const & p4 );
template< typename T > void dihedral_degrees( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & p3, xyzVector< T > const & p4, T & angle );
template< typename T > T dihedral_degrees( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & p3, xyzVector< T > const & p4 );
template< typename T > void dihedral( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & p3, xyzVector< T > const & p4, T & angle );
template< typename T > T dihedral( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & p3, xyzVector< T > const & p4 );
template< typename T > xyzMatrix< T > rotation_matrix( xyzVector< T > const & axis, T const & theta );
template< typename T > xyzVector< T > rotation_axis( xyzMatrix< T > const & R, T & theta );
template< typename T > xyzVector< T > eigenvalue_jacobi( xyzMatrix< T > const & a, T const & tol );
template< typename T > xyzVector< T > eigenvector_jacobi( xyzMatrix< T > const & a, T const & tol, xyzMatrix< T > & J );
template< typename T > void jacobi_rotation( xyzMatrix< T > const & m, int const i, int const j, xyzMatrix< T > & r );
template< typename T > sphericalVector< T > xyz_to_spherical(xyzVector< T > const & xyz );
template< typename T > xyzVector< T > spherical_to_xyz(sphericalVector< T > const & spherical);


} // namespace numeric


#endif // INCLUDED_numeric_xyz_functions_FWD_HH

