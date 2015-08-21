// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyz.functions.hh
/// @brief  xyzVector and xyzMatrix functions
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_xyz_functions_hh
#define INCLUDED_numeric_xyz_functions_hh


// Package headers
#include <numeric/conversions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/sphericalVector.hh>


//utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cmath>
#include <cstdlib>
#include <vector>

//ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

namespace numeric {


// Products


/// @brief xyzMatrix * xyzVector
/// @note  Same as product( xyzMatrix, xyzVector )
//template< typename T >
//inline
//xyzVector< T >
//operator *( xyzMatrix< T > const & m, xyzVector< T > const & v )
//{
// return xyzVector< T >(
//  m.xx_ * v.x_ + m.xy_ * v.y_ + m.xz_ * v.z_,
//  m.yx_ * v.x_ + m.yy_ * v.y_ + m.yz_ * v.z_,
//  m.zx_ * v.x_ + m.zy_ * v.y_ + m.zz_ * v.z_
// );
//}

/// @brief return the point closest to point p3 that lies on the line
/// defined by p1 and p2
template< typename T >
inline
xyzVector< T >
closest_point_on_line( xyzVector< T > const & p1, xyzVector< T > const & p2, xyzVector< T > const & q )
{
	xyzVector<T> u = p2-p1;
	xyzVector<T> pq = q-p1;
	xyzVector<T> w2 = pq-(u*(dot_product(pq,u)/u.magnitude_squared()));
	xyzVector<T> point = q-w2;
	return point;
}

/// @brief calculate center of mass for coordinates
template< typename T >
inline
xyzVector< T >
center_of_mass( utility::vector1< xyzVector< T > > const & coords )
{
	xyzVector< T > center_of_mass( 0.0, 0.0, 0.0 );
	for ( typename utility::vector1< xyzVector< T > >::const_iterator it = coords.begin();
			it != coords.end();
			++it ) {
		center_of_mass += *it;
	}
	center_of_mass /= coords.size();
	return center_of_mass;
}

/// @brief xyzMatrix * xyzVector product
/// @note  Same as xyzMatrix * xyzVector
template< typename T >
inline
xyzVector< T >
product( xyzMatrix< T > const & m, xyzVector< T > const & v )
{
	return xyzVector< T >(
		m.xx_ * v.x_ + m.xy_ * v.y_ + m.xz_ * v.z_,
		m.yx_ * v.x_ + m.yy_ * v.y_ + m.yz_ * v.z_,
		m.zx_ * v.x_ + m.zy_ * v.y_ + m.zz_ * v.z_
	);
}


/// @brief xyzMatrix * xyzVector in-place product
/// @note  Input xyzVector is modified
template< typename T >
inline
xyzVector< T > &
inplace_product( xyzMatrix< T > const & m, xyzVector< T > & v )
{
	T const x = m.xx_ * v.x_ + m.xy_ * v.y_ + m.xz_ * v.z_;
	T const y = m.yx_ * v.x_ + m.yy_ * v.y_ + m.yz_ * v.z_;
	v.z_      = m.zx_ * v.x_ + m.zy_ * v.y_ + m.zz_ * v.z_;
	v.x_ = x;
	v.y_ = y;
	return v;
}


/// @brief xyzMatrix^T * xyzVector product
template< typename T >
inline
xyzVector< T >
transpose_product( xyzMatrix< T > const & m, xyzVector< T > const & v )
{
	return xyzVector< T >(
		m.xx_ * v.x_ + m.yx_ * v.y_ + m.zx_ * v.z_,
		m.xy_ * v.x_ + m.yy_ * v.y_ + m.zy_ * v.z_,
		m.xz_ * v.x_ + m.yz_ * v.y_ + m.zz_ * v.z_
	);
}


/// @brief xyzMatrix^T * xyzVector in-place transpose product
/// @note  Input xyzVector is modified
template< typename T >
inline
xyzVector< T > &
inplace_transpose_product( xyzMatrix< T > const & m, xyzVector< T > & v )
{
	T const x = m.xx_ * v.x_ + m.yx_ * v.y_ + m.zx_ * v.z_;
	T const y = m.xy_ * v.x_ + m.yy_ * v.y_ + m.zy_ * v.z_;
	v.z_      = m.xz_ * v.x_ + m.yz_ * v.y_ + m.zz_ * v.z_;
	v.x_ = x;
	v.y_ = y;
	return v;
}


/// @brief xyzVector xyzVector outer product
template< typename T >
inline
xyzMatrix< T >
outer_product( xyzVector< T > const & a, xyzVector< T > const & b )
{
	return xyzMatrix< T >(
		a.x_ * b.x_, a.x_ * b.y_, a.x_ * b.z_,
		a.y_ * b.x_, a.y_ * b.y_, a.y_ * b.z_,
		a.z_ * b.x_, a.z_ * b.y_, a.z_ * b.z_
	);
}

//Averages


/// @brief geometric center
/// @note compute the geometric center of a list of points
/*
template<typename T>
inline
xyzVector<T>
geometric_center(std::vector<xyzVector<T> > const & points)
{
xyzVector<T> total(0.0,0.0,0.0);
std::vector<xyzVector<T> >::iterator points_it = points.begin();
for(; points_it != points.end();++points_it)
{
total += *points_it;
}
return total/points.size();
}
*/

// Projections


/// @brief Projection matrix onto the line through a vector
template< typename T >
inline
xyzMatrix< T >
projection_matrix( xyzVector< T > const & v )
{
	return ( xyzMatrix< T >(
		v.x_ * v.x_, v.x_ * v.y_, v.x_ * v.z_,
		v.y_ * v.x_, v.y_ * v.y_, v.y_ * v.z_,
		v.z_ * v.x_, v.z_ * v.y_, v.z_ * v.z_
		) /= v.length_squared() );
}

template< typename T >
inline
xyzMatrix< T >
inverse( xyzMatrix< T > const & a ) {
	T D = a.det();
	return xyzMatrix< T >(
		(a.yy_*a.zz_-a.yz_*a.zy_)/D, -(a.xy_*a.zz_-a.xz_*a.zy_)/D,  (a.xy_*a.yz_-a.xz_*a.yy_)/D,
		-(a.yx_*a.zz_-a.zx_*a.yz_)/D,  (a.xx_*a.zz_-a.xz_*a.zx_)/D, -(a.xx_*a.yz_-a.xz_*a.yx_)/D,
		(a.yx_*a.zy_-a.zx_*a.yy_)/D, -(a.xx_*a.zy_-a.xy_*a.zx_)/D,  (a.xx_*a.yy_-a.xy_*a.yx_)/D
	);
}


// Angles


/// @brief Plane angle in radians: angle value passed
/// @note  Given thre positions in a chain ( p1, p2, p3 ), calculates the plane
///   angle in radians between the vectors p2->p1 and p2->p3
/// @note  Angle returned is on [ 0, pi ]
template< typename T >
inline
void
angle_radians(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	T & angle // Angle (radians)
)
{
	// Two normalized directional vectors hold the relative directions of
	// p1 and p3 relative to p2
	xyzVector< T > const a( ( p2 - p1 ).normalize_or_zero() );
	xyzVector< T > const b( ( p2 - p3 ).normalize_or_zero() );

	angle = std::acos( sin_cos_range( dot(a, b) ) );
}


/// @brief Plane angle in radians: angle value returned
/// @note  Given three positions in a chain ( p1, p2, p3 ), calculates the plane
///   angle in radians between the vectors p2->p1 and p2->p3
/// @note  Angle returned is on [ 0, pi ]
template< typename T >
inline
T // Angle (radians)
angle_radians(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3
)
{
	T angle;
	angle_radians( p1, p2, p3, angle );
	return angle;
}

// PyRosetta, creating concreate function for template one
inline double angle_radians_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3) { return angle_radians(p1, p2, p3); }


/// @brief Plane angle in degrees: angle value returned
/// @note  Given three positions in a chain ( p1, p2, p3 ), calculates the plane
///   angle in degrees between the vectors p2->p1 and p2->p3
/// @note  Angle returned is on [ 0, 180 ]
template< typename T >
inline
T // Angle (degrees)
angle_degrees(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3
)
{
	T angle;
	angle_radians( p1, p2, p3, angle );
	return conversions::degrees(angle);
}

// PyRosetta, creating concreate function for template one
inline double angle_degrees_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3) { return angle_degrees(p1, p2, p3); }


/// @brief Angle between two vectors in radians
/// @note  Given two vectors (p1->p2 & p3->p4),
/// calculate the angle between them
/// @note  Angle returned is on [ 0, pi ]
template< typename T >
inline
T // Angle (radians)
angle_radians(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4
)
{
	xyzVector< T > const a( ( p2 - p1 ).normalize_or_zero() );
	xyzVector< T > const b( ( p4 - p3 ).normalize_or_zero() );

	T angle = std::acos( sin_cos_range( dot(a, b) ) );
	return angle;
}

// PyRosetta, creating concreate function for template one
inline double angle_radians_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4) { return angle_radians(p1, p2, p3, p4); }


/// @brief Angle between two vectors in radians
/// @note  Given two vectors (p1->p2 & p3->p4),
/// calculate the angle between them
/// @note  Angle returned is on [ 0, pi ]
template< typename T >
inline
T // Angle (radians)
angle_degrees(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4
)
{
	T angle = angle_radians( p1, p2, p3, p4 );
	return conversions::degrees(angle);
}

// PyRosetta, creating concreate function for template one
inline double angle_degrees_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4) { return angle_degrees(p1, p2, p3, p4); }


/// @brief Dihedral (torsion) angle in radians: angle value passed
/// @note  Given four positions in a chain ( p1, p2, p3, p4 ), calculates the dihedral
///   (torsion) angle in radians between the vectors p2->p1 and p3->p4 while sighting
///   along the axis defined by the vector p2->p3 (positive indicates right handed twist)
/// @note  Angle returned is on [ -pi, pi ]
/// @note  Degenerate cases are handled and assigned a zero angle but assumed rare
///   (wrt performance tuning)
/// @note For a reference on the determination of the dihedral angle formula see:
///   http://www.math.fsu.edu/~quine/IntroMathBio_04/torsion_pdb/torsion_pdb.pdf
template< typename T >
inline
void
dihedral_radians(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4,
	T & angle // Angle (radians)
)
{
	static T const ZERO( 0 );

	// Three normalized directional vectors hold the relative directions of
	// the points (the dihedral angle formula used below is only valid for
	// normalized vectors)
	xyzVector< T > const a( ( p2 - p1 ).normalize_or_zero() );
	xyzVector< T > const b( ( p3 - p2 ).normalize_or_zero() );
	xyzVector< T > const c( ( p4 - p3 ).normalize_or_zero() );

	// Compute the dihedral angle: Degenerate cases are assigned a zero angle by definition
	// Degenerate cases: Coincident adjacent points or parallel adjacent directional vectors
	if ( ( b == ZERO ) ) { // p2 == p3: Handle specially: atan2( 0, -dot( a, c ) ) == 0 or pi
		angle = ZERO;
	} else { // Use the formula
		// Degenerate cases a == 0, c == 0, a||b, and b||c give x == y == 0
		T const x = -dot( a, c ) + ( dot( a, b ) * dot( b, c ) );
		T const y = dot( a, cross( b, c ) );
		// Angle in [ -pi, pi ]
		angle = ( ( y != ZERO ) || ( x != ZERO ) ? std::atan2( y, x ) : ZERO );
	}
}

// PyRosetta, creating concreate function for template one
inline void dihedral_radians_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4,
	double & angle // Angle (radians)
) { dihedral_radians(p1, p2, p3, p4, angle); }


/// @brief Dihedral (torsion) angle in radians: angle value returned
/// @note  Given four positions in a chain ( p1, p2, p3, p4 ), calculates the dihedral
///   (torsion) angle in radians between the vectors p2->p1 and p3->p4 while sighting
///   along the axis defined by the vector p2->p3 (positive indicates right handed twist)
/// @note  Angle returned is on [ -pi, pi ]
/// @note  Degenerate cases are handled and assigned a zero angle but assumed rare
///   (wrt performance tuning)
/// @note For a reference on the determination of the dihedral angle formula see:
///   http://www.math.fsu.edu/~quine/IntroMathBio_04/torsion_pdb/torsion_pdb.pdf
template< typename T >
inline
T // Angle (radians)
dihedral_radians(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4
)
{
	T angle;
	dihedral_radians( p1, p2, p3, p4, angle );
	return angle;
}

// PyRosetta, creating concreate function for template one
inline double dihedral_radians_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4
) { return dihedral_radians(p1, p2, p3, p4); }


/// @brief Dihedral (torsion) angle in degrees: angle value passed
template< typename T >
inline
void
dihedral_degrees(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4,
	T & angle // Angle (degrees)
)
{
	dihedral_radians( p1, p2, p3, p4, angle );
	conversions::to_degrees( angle );
}

// PyRosetta, creating concreate function for template one
inline void dihedral_degrees_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4,
	double & angle // Angle (radians)
) { dihedral_degrees(p1, p2, p3, p4, angle); }


/// @brief Dihedral (torsion) angle in degrees: angle value returned
template< typename T >
inline
T // Angle (degrees)
dihedral_degrees(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4
)
{
	return conversions::degrees( dihedral_radians( p1, p2, p3, p4 ) );
}


// PyRosetta, creating concreate function for template one
inline double dihedral_degrees_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4
) { return dihedral_degrees(p1, p2, p3, p4); }


/// @brief Dihedral (torsion) angle in degrees: angle value passed
/// @note  This is a Rosetta++ compatibility version that operates in degrees
template< typename T >
inline
void
dihedral(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4,
	T & angle // Angle (degrees)
)
{
	dihedral_radians( p1, p2, p3, p4, angle );
	conversions::to_degrees( angle );
}

// PyRosetta, creating concreate function for template one
inline void dihedral_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4,
	double & angle // Angle (degrees)
) { dihedral(p1, p2, p3, p4, angle); }


/// @brief Dihedral (torsion) angle in degrees: angle value returned
/// @note  This is a Rosetta++ compatibility version that operates in degrees
template< typename T >
inline
T // Angle (degrees)
dihedral(
	xyzVector< T > const & p1,
	xyzVector< T > const & p2,
	xyzVector< T > const & p3,
	xyzVector< T > const & p4
)
{
	return conversions::degrees( dihedral_radians( p1, p2, p3, p4 ) );
}


// PyRosetta, creating concreate function for template one
inline double dihedral_double(
	xyzVector< double > const & p1,
	xyzVector< double > const & p2,
	xyzVector< double > const & p3,
	xyzVector< double > const & p4
) { return dihedral(p1, p2, p3, p4); }


// @brief Rotation matrix for rotation from axis-angle representation
// @note Magnitude of rotation (in radians) is taken as axis_angle.magnitude
template< typename T >
inline
xyzMatrix< T >
rotation_matrix(xyzVector< T > const & axis_angle)
{
	// rotation_matrix performs axis vector normalization
	return rotation_matrix(axis_angle, axis_angle.magnitude());
}

/// @brief Rotation matrix for rotation about an axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
rotation_matrix(
	xyzVector< T > const & axis,
	T const & theta // Angle (radians)
)
{
	xyzVector< T > const n( axis.normalized() );

	T const sin_theta( std::sin( theta ) );
	T const cos_theta( std::cos( theta ) );

	xyzMatrix< T > R( projection_matrix( n ) *= ( T( 1 ) - cos_theta ) );

	R.xx_ += cos_theta;        R.xy_ -= sin_theta * n.z_; R.xz_ += sin_theta * n.y_;
	R.yx_ += sin_theta * n.z_; R.yy_ += cos_theta;        R.yz_ -= sin_theta * n.x_;
	R.zx_ -= sin_theta * n.y_; R.zy_ += sin_theta * n.x_; R.zz_ += cos_theta;

	return R;
}


/// @brief Rotation matrix for rotation about an axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
rotation_matrix_radians(
	xyzVector< T > const & axis,
	T const & theta // Angle (radians)
)
{
	return rotation_matrix( axis, theta );
}


/// @brief Rotation matrix for rotation about an axis by an angle in degrees
template< typename T >
inline
xyzMatrix< T >
rotation_matrix_degrees(
	xyzVector< T > const & axis,
	T const & theta // Angle (degrees)
)
{
	return rotation_matrix( axis, conversions::radians( theta ) );
}


/// @brief Rotation matrix for rotation about the x axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
x_rotation_matrix(
	T const & theta // Angle (radians)
)
{
	T const sin_theta( std::sin( theta ) );
	T const cos_theta( std::cos( theta ) );

	return xyzMatrix< T >::rows(
		T( 1 ), T( 0 ),     T( 0 ),
		T( 0 ), cos_theta, -sin_theta,
		T( 0 ), sin_theta,  cos_theta
	);
}


/// @brief Rotation matrix for rotation about the x axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
x_rotation_matrix_radians(
	T const & theta // Angle (radians)
)
{
	return x_rotation_matrix( theta );
}


/// @brief Rotation matrix for rotation about the x axis by an angle in degrees
template< typename T >
inline
xyzMatrix< T >
x_rotation_matrix_degrees(
	T const & theta // Angle (degrees)
)
{
	return x_rotation_matrix( conversions::radians( theta ) );
}


/// @brief Rotation matrix for rotation about the y axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
y_rotation_matrix(
	T const & theta // Angle (radians)
)
{
	T const sin_theta( std::sin( theta ) );
	T const cos_theta( std::cos( theta ) );

	return xyzMatrix< T >::rows(
		cos_theta, T( 0 ), sin_theta,
		T( 0 ),    T( 1 ), T( 0 ),
		-sin_theta, T( 0 ), cos_theta
	);
}


/// @brief Rotation matrix for rotation about the y axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
y_rotation_matrix_radians(
	T const & theta // Angle (radians)
)
{
	return y_rotation_matrix( theta );
}


/// @brief Rotation matrix for rotation about the y axis by an angle in degrees
template< typename T >
inline
xyzMatrix< T >
y_rotation_matrix_degrees(
	T const & theta // Angle (degrees)
)
{
	return y_rotation_matrix( conversions::radians( theta ) );
}


/// @brief Rotation matrix for rotation about the z axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
z_rotation_matrix(
	T const & theta // Angle (radians)
)
{
	T const sin_theta = static_cast< T > ( std::sin( theta ) );
	T const cos_theta = static_cast< T > ( std::cos( theta ) );

	return xyzMatrix< T >::rows(
		cos_theta, -sin_theta, T( 0 ),
		sin_theta,  cos_theta, T( 0 ),
		T( 0 )   ,  T( 0 ),    T( 1 )
	);
}


/// @brief Rotation matrix for rotation about the z axis by an angle in radians
template< typename T >
inline
xyzMatrix< T >
z_rotation_matrix_radians(
	T const & theta // Angle (radians)
)
{
	return z_rotation_matrix( theta );
}


/// @brief Rotation matrix for rotation about the z axis by an angle in degrees
template< typename T >
inline
xyzMatrix< T >
z_rotation_matrix_degrees(
	T const & theta // Angle (degrees)
)
{
	return z_rotation_matrix( conversions::radians( theta ) );
}


/// @brief  Helper function to find the rotation to optimally transform the vectors A1-B1 to vectors A2-B2
template< typename T >
inline
xyzMatrix< T >
alignVectorSets(
	xyzVector< T > A1,
	xyzVector< T > B1,
	xyzVector< T > A2,
	xyzVector< T > B2 ) {
	// 1) find rotation to align canonic +z to each vector pair's +z (defined as the avg of the two vectors)
	xyzVector< T > X1 = (A1+B1); X1.normalize();
	xyzVector< T > X2 = (A2+B2); X2.normalize();

	T cb1 = X1[2], sb1 = sqrt(1-(X1[2]*X1[2]));
	T Rg1 = sqrt((X1[0]*X1[0])+(X1[1]*X1[1]));
	T cg1 = X1[1]/Rg1, sg1 = X1[0]/Rg1;

	T cb2 = X2[2], sb2 = sqrt(1-(X2[2]*X2[2]));
	T Rg2 = sqrt((X2[0]*X2[0])+(X2[1]*X2[1]));
	T cg2 = X2[1]/Rg2, sg2 = X2[0]/Rg2;

	xyzMatrix< T > R1gb, R2gb, R1gb_i, R2gb_i;
	R1gb.xx( cg1 ); R1gb.xy(cb1*sg1); R1gb.xz(sb1*sg1);
	R1gb.yx(-sg1 ); R1gb.yy(cb1*cg1); R1gb.yz(sb1*cg1);
	R1gb.zx( 0   ); R1gb.zy(-sb1)   ; R1gb.zz(cb1);

	R2gb.xx( cg2 ); R2gb.xy(cb2*sg2); R2gb.xz(sb2*sg2);
	R2gb.yx(-sg2 ); R2gb.yy(cb2*cg2); R2gb.yz(sb2*cg2);
	R2gb.zx( 0   ); R2gb.zy(-sb2)   ; R2gb.zz(cb2);

	R1gb_i = numeric::inverse( R1gb );
	R2gb_i = numeric::inverse( R2gb );

	// 2) now choose one of the two vectors (A1/A2) to define the xz plane
	A1.normalize();
	A2.normalize();
	numeric::xyzVector< T > RgbA1 = R1gb_i*A1;
	numeric::xyzVector< T > RgbA2 = R2gb_i*A2;

	T Ra1 = sqrt((RgbA1[0]*RgbA1[0])+(RgbA1[1]*RgbA1[1]));
	T Ra2 = sqrt((RgbA2[0]*RgbA2[0])+(RgbA2[1]*RgbA2[1]));
	T ca1 = RgbA1[1]/Ra1, sa1 = RgbA1[0]/Ra1;
	T ca2 = RgbA2[1]/Ra2, sa2 = RgbA2[0]/Ra2;

	numeric::xyzMatrix< T > R1, R2, R1_i, R;
	R1.xx( -sa1*cb1*sg1 + ca1*cg1 ); R1.xy(  ca1*cb1*sg1 + sa1*cg1 ); R1.xz(  sb1*sg1 );
	R1.yx( -sa1*cb1*cg1 - ca1*sg1 ); R1.yy(  ca1*cb1*cg1 - sa1*sg1 ); R1.yz(  sb1*cg1 );
	R1.zx(  sa1*sb1 );               R1.zy( -ca1*sb1 );               R1.zz(  cb1 );

	R2.xx( -sa2*cb2*sg2 + ca2*cg2 ); R2.xy(  ca2*cb2*sg2 + sa2*cg2 ); R2.xz(  sb2*sg2 );
	R2.yx( -sa2*cb2*cg2 - ca2*sg2 ); R2.yy(  ca2*cb2*cg2 - sa2*sg2 ); R2.yz(  sb2*cg2 );
	R2.zx(  sa2*sb2 );               R2.zy( -ca2*sb2 );               R2.zz(  cb2 );

	// 3) the rotation matrix first rotates 1 to canonical, then canonical to 2
	R1_i = numeric::inverse( R1 );
	R = R2*R1_i;

	return R;
}

/// @brief Transformation from rotation matrix to magnitude of helical rotation
/// @note  Input matrix must be orthogonal
/// @note  Orientation of axis chosen so that the angle of rotation is non-negative [0,pi]
//  @note  numeric::rotation_axis returns both axis and angle of rotation
template< typename T >
inline
T
rotation_angle( xyzMatrix< T > const & R)
{
	using std::abs;
	using std::acos;
	using std::sqrt;

	// This would be good here but slow and unclear what tolerance to use
	// assert( rm.orthogonal() );

	static T const ZERO( 0 );
	static T const ONE( 1 );
	static T const TWO( 2 );

	T const cos_theta = sin_cos_range( ( R.trace() - ONE ) / TWO );

	T const tolerance = NumericTraits< T >::sin_cos_tolerance();
	if ( cos_theta > -ONE + tolerance && cos_theta < ONE - tolerance ) {
		return acos( cos_theta );
	} else if ( cos_theta >= ONE - tolerance ) {
		// R is the identity matrix, return an arbitrary axis of rotation
		return ZERO;
	} else { // cos_theta <= -ONE + tolerance
		return NumericTraits< T >::pi(); // theta == pi
	}
}

/// @brief Transformation from rotation matrix to helical axis of rotation
/// @note  Input matrix must be orthogonal
/// @note  Angle of rotation is also returned
/// @note  Orientation of axis chosen so that the angle of rotation is non-negative [0,pi]
template< typename T >
inline
xyzVector< T >
rotation_axis( xyzMatrix< T > const & R, T & theta )
{
	using std::abs;
	using std::acos;
	using std::sqrt;

	// This would be good here but slow and unclear what tolerance to use
	// assert( rm.orthogonal() );

	static T const ZERO( 0 );
	static T const ONE( 1 );
	static T const TWO( 2 );

	T const cos_theta = sin_cos_range( ( R.trace() - ONE ) / TWO );

	T const tolerance = NumericTraits< T >::sin_cos_tolerance();
	if ( cos_theta > -ONE + tolerance && cos_theta < ONE - tolerance ) {
		// Compute sign and absolute value of axis vector elements from matrix elements
		// Sign of axis vector is chosen to correspond to a positive sin_theta value
		T x = ( R.zy_ > R.yz_ ? ONE : -ONE ) *
			sqrt( max( ZERO, ( R.xx_ - cos_theta ) / ( ONE - cos_theta ) ) );
		T y = ( R.xz_ > R.zx_ ? ONE : -ONE ) *
			sqrt( max( ZERO, ( R.yy_ - cos_theta ) / ( ONE - cos_theta ) ) );
		T z = ( R.yx_ > R.xy_ ? ONE : -ONE ) *
			sqrt( max( ZERO, ( R.zz_ - cos_theta ) / ( ONE - cos_theta ) ) );
		// Above method appears to cover a greater range of cases than the original method:
		//
		// return ( xyzVector< T >( R.zy_ - R.yz_, R.xz_ - R.zx_, R.yx_ - R.xy_ )
		//          /= ( TWO * sine_theta ) );
		//
		// and is more stable for small sine_theta

		theta = acos( cos_theta );
		assert( std::abs( x*x + y*y + z*z - 1 ) <= T( 0.01 ) );

		return xyzVector< T >( x, y, z );
	} else if ( cos_theta >= ONE - tolerance ) {
		// R is the identity matrix, return an arbitrary axis of rotation
		theta = ZERO;
		return xyzVector< T >( ONE, ZERO, ZERO );
	} else { // cos_theta <= -ONE + tolerance
		// R is of the form 2*n*n^T - I, theta == pi
		xyzMatrix< T > const nnT( xyzMatrix< T >( R ).add_diagonal( ONE, ONE, ONE ) /= TWO );
		T x, y, z;
		// Since theta = pi, both n and -n are equally valid axis vectors
		// Use convention to take first non-zero coordinate positive
		if ( nnT.xx_ > ZERO + tolerance ) {
			// Note: x is positive, but negative would also work
			x = sqrt( nnT.xx_ );
			y = nnT.yx_ / x;
			z = nnT.zx_ / x;
		} else if ( nnT.yy_ > ZERO + tolerance ) {
			// nnT.xx_ = n.x_ == ZERO, in this case
			x = ZERO;
			// Note: y is positive, but negative would also work
			y = sqrt( nnT.yy_ );
			z = nnT.zy_ / y;
		} else { // nnT.yy_ = n.y_ == ZERO, also, in this case
			// If we get here nnT.zz_ must be positive!
			assert( nnT.zz_ > ZERO + tolerance );
			x = ZERO;
			y = ZERO;
			// Note: z is positive, but negative would also work
			z = sqrt( nnT.zz_ );
		}
		theta = NumericTraits< T >::pi(); // theta == pi
		// For a valid orthogonal matrix R, axis should be a normal vector
		assert( std::abs( x*x + y*y + z*z - 1 ) <= T( 0.01 ) );
		return xyzVector< T >( x, y, z );
	}
}

/// @brief Transformation from rotation matrix to compact axis-angle representation
/// @note  Input matrix must be orthogonal
/// @note  Orientation of axis chosen so that the angle of rotation is non-negative [0,pi]
//  @note  Resulting vector will be oriented in axis of rotation with magnitude equal to magnitude of rotation.
template< typename T>
inline
xyzVector< T >
rotation_axis_angle(xyzMatrix< T > const & R)
{
	T theta;
	xyzVector< T > vec = rotation_axis(R, theta);

	return vec * theta;
}

// Eigenvalues/vectors


/// @brief Classic Jacobi algorithm for the eigenvalues of a real symmetric matrix
/// @note  Use eigenvector_jacobi if eigenvectors are also desired
template< typename T >
inline
xyzVector< T >
eigenvalue_jacobi( xyzMatrix< T > const & a, T const & tol )
{
	using std::abs;

	// May need a tolerance based test here
	assert( ( a.xy_ == a.yx_ ) && ( a.xz_ == a.zx_ ) && ( a.yz_ == a.zy_ ) );

	// Copy matrix as it will be modified by the algorithm
	xyzMatrix< T > m( a );

	// Initialize the off-diagonal, upper-triangular sum of squares
	T off = ( m.xy_ * m.xy_ ) + ( m.xz_ * m.xz_ ) + ( m.yz_ * m.yz_ );

	int i, j, n_iterations = 0;
	xyzMatrix< T > r;
	while ( off > tol ) {
		++n_iterations;

		// Ensure number of iterations does not exceed 50:
		// otherwise, re-evaluate stopping criterion
		assert( n_iterations <= 50 );

		// Determine index of upper-triangular element that will be zeroed out
		if ( std::abs( m.xy_ ) >= std::abs( m.xz_ ) ) {
			if ( std::abs( m.xy_ ) >= std::abs( m.yz_ ) ) {
				i = 1; j = 2;
			} else {
				i = 2; j = 3;
			}
		} else if ( std::abs( m.xz_ ) >= std::abs( m.yz_ ) ) {
			i = 1; j = 3;
		} else {
			i = 2; j = 3;
		}

		// After four iterations, skip the rotation if the off-diagonal element is small
		T const ij_scaled = std::abs( T( 100 ) * m(i,j) );
		if ( ( n_iterations > 4 )
				&& std::abs( m(i,i) ) + ij_scaled == std::abs( m(i,i) )
				&& std::abs( m(j,j) ) + ij_scaled == std::abs( m(j,j) ) ) {
			m(i,j) = m(j,i) = T( 0 );
		} else {
			// Compute the rotation matrix
			jacobi_rotation( m, i, j, r );

			// Zero out the i,j and j,i elements
			m.right_multiply_by( r );
			m.left_multiply_by_transpose( r );
		}

		// Recalculate the off-diagonal, upper-triangular sum of squares
		off = ( m.xy_ * m.xy_ ) + ( m.xz_ * m.xz_ ) + ( m.yz_ * m.yz_ );
	}

	return xyzVector< T >( m.xx_, m.yy_, m.zz_ );
}


/// @brief Classic Jacobi algorithm for the eigenvalues and eigenvectors of a
///   real symmetric matrix
/// @note  Use eigenvalue_jacobi if eigenvectors are not desired
template< typename T >
inline
xyzVector< T >
eigenvector_jacobi( xyzMatrix< T > const & a, T const & tol, xyzMatrix< T > & J )
{
	using std::abs;

	// May need a tolerance based test here
	assert( ( a.xy_ == a.yx_ ) && ( a.xz_ == a.zx_ ) && ( a.yz_ == a.zy_ ) );

	// Copy matrix as it will be modified by the algorithm
	xyzMatrix< T > m( a );

	// Initialize the off-diagonal, upper-triangular sum of squares
	T off = ( m.xy_ * m.xy_ ) + ( m.xz_ * m.xz_ ) + ( m.yz_ * m.yz_ );

	J.to_identity();
	int i, j, n_iterations = 0;
	xyzMatrix< T > r;
	while ( off > tol ) {
		++n_iterations;

		// Ensure number of iterations does not exceed 50:
		// otherwise, re-evaluate stopping criterion
		assert( n_iterations <= 50 );

		// Determine index of upper-triangular element that will be zeroed out
		if ( std::abs( m.xy_ ) >= std::abs( m.xz_ ) ) {
			if ( std::abs( m.xy_ ) >= std::abs( m.yz_ ) ) {
				i = 1; j = 2;
			} else {
				i = 2; j = 3;
			}
		} else if ( std::abs( m.xz_ ) >= std::abs( m.yz_ ) ) {
			i = 1; j = 3;
		} else {
			i = 2; j = 3;
		}

		// After four iterations, skip the rotation if the off-diagonal element is small
		T const ij_scaled = std::abs( T( 100 ) * m(i,j) );
		if ( ( n_iterations > 4 )
				&& std::abs( m(i,i) ) + ij_scaled == std::abs( m(i,i) )
				&& std::abs( m(j,j) ) + ij_scaled == std::abs( m(j,j) ) ) {
			m(i,j) = m(j,i) = T( 0 );
		} else {
			// Compute the rotation matrix
			jacobi_rotation( m, i, j, r );

			// Zero out the i,j and j,i elements
			m.right_multiply_by( r );
			m.left_multiply_by_transpose( r );

			// Accumulate the rotation transformations to form the matrix of eigenvectors
			J *= r;
		}

		// Recalculate the off-diagonal, upper-triangular sum of squares
		off = ( m.xy_ * m.xy_ ) + ( m.xz_ * m.xz_ ) + ( m.yz_ * m.yz_ );
	}

	return xyzVector< T >( m.xx_, m.yy_, m.zz_ );
}


/// @brief Jacobi rotation
/// @note  Compute the orthogonal transformation used to zero out a pair of
///        off-diagonal elements
template< typename T >
inline
void
jacobi_rotation( xyzMatrix< T > const & m, int const i, int const j, xyzMatrix< T > & r )
{
	using std::abs;
	using std::sqrt;

	assert( ( i > 0 ) && ( i <= 3 ) );
	assert( ( j > 0 ) && ( j <= 3 ) );
	assert( i != j );

	static T const ZERO( 0 );
	static T const ONE( 1 );

	T const tau = ( m(i,i) - m(j,j) ) / ( 2 * m(i,j) );
	T const t = ( tau < ZERO ? -ONE : ONE ) / ( std::abs( tau ) + sqrt( ONE + ( tau * tau ) ) );

	T const c = ONE / sqrt( ONE + ( t * t ) );
	T const s = c * t;

	r.to_identity();
	r(i,i) = c; r(i,j) = -s;
	r(j,i) = s; r(j,j) = c;
}

template<typename T>
inline
sphericalVector<T>
xyz_to_spherical(xyzVector<T> const & xyz)
{

	sphericalVector<T> spherical;
	spherical.radius(sqrt((xyz.x()*xyz.x())+(xyz.y()*xyz.y())+(xyz.z()*xyz.z())));
	spherical.theta(acos(xyz.z()/spherical.radius())*(1/numeric::constants::f::pi_over_180)) ;
	spherical.phi(atan2((xyz.y()),(xyz.x()))*(1/numeric::constants::f::pi_over_180));
	return spherical;
}

template<typename T>
inline
xyzVector<T>
spherical_to_xyz(sphericalVector<T> const & spherical)
{
	xyzVector<T> xyz;
	xyz.x(spherical.radius()*sin(spherical.theta()*numeric::constants::f::pi_over_180)*cos(spherical.phi()*numeric::constants::f::pi_over_180));
	xyz.y(spherical.radius()*sin(spherical.theta()*numeric::constants::f::pi_over_180)*sin(spherical.phi()*numeric::constants::f::pi_over_180));
	xyz.z(spherical.radius()*cos(spherical.theta()*numeric::constants::f::pi_over_180));
	return xyz;
}


/// @brief convert a string of comma separated values "0.2,0.4,0.3" to an xyzVector
template<typename T>
inline
xyzVector<T>
comma_seperated_string_to_xyz(std::string triplet)
{
	utility::vector1<std::string> split_string(utility::string_split(triplet, ','));
	runtime_assert(split_string.size() == 3);
	xyzVector<T> xyz;
	xyz.x(utility::from_string(split_string[1],T(0.0)));
	xyz.y(utility::from_string(split_string[2],T(0.0)));
	xyz.z(utility::from_string(split_string[3],T(0.0)));
	return xyz;

}


/// @brief convert a vector1 of xyzVectors to an FArray2D
template<typename T>
inline
ObjexxFCL::FArray2D<T>
vector_of_xyzvectors_to_FArray(utility::vector1< xyzVector<T> > const & input)
{
	ObjexxFCL::FArray2D< T > output(3,input.size());
	for ( numeric::Real index = 1; index <= input.size(); ++index ) {
		output(1,(int)index) = input[(int)index].x(); // bazzoli: added cast to silence warning.
		output(2,(int)index) = input[(int)index].y();
		output(3,(int)index) = input[(int)index].z();
	}
	return output;
}

/// @brief convert an FArray2D to a vector of xyzVectors
template<typename T>
inline
utility::vector1< xyzVector<T> >
FArray_to_vector_of_xyzvectors(ObjexxFCL::FArray2D<T> const & input)
{
	assert(input.size1() == 3);
	utility::vector1< xyzVector<T> > output(input.size2(),xyzVector<T>());
	for ( numeric::Real index = 1; index <= input.size2(); ++index ) {
		output[(int)index].x(input(1,(int)index));
		output[(int)index].y(input(2,(int)index));
		output[(int)index].z(input(3,(int)index));
	}

	return output;
}

/// @brief convert a 3x3 FArray 2D to an xyzMatrix
template<typename T>
inline
numeric::xyzMatrix<T>
FArray_to_xyzmatrix(ObjexxFCL::FArray2D<T> const & input)
{
	assert(input.size1() == 3 && input.size2() == 3);

	return xyzMatrix<T>::rows(
		input(1,1),input(1,2),input(1,3),
		input(2,1),input(2,2),input(2,3),
		input(3,1),input(3,2),input(3,3)
	);
}

/// @brief convert an xyzMatrix to a 3x3 FArray 2D
template<typename T>
inline
ObjexxFCL::FArray2D<T> xyzmatrix_to_FArray(numeric::xyzMatrix<T> const & input)
{
	ObjexxFCL::FArray2D<T> output(3,3);
	output(1,1) = input.xx();
	output(1,2) = input.xy();
	output(1,3) = input.xz();
	output(2,1) = input.yx();
	output(2,2) = input.yy();
	output(2,3) = input.yz();
	output(3,1) = input.zx();
	output(3,2) = input.zy();
	output(3,3) = input.zz();

	return output;

}

} // namespace numeric


#endif // INCLUDED_numeric_xyz_functions_HH
