// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/HomogeneousTransform.hh
/// @brief  Fast coordinate frame container
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
///
/// @remarks
///  @li Inline, loop-free functions for speed
///  @li Non-virtual destructor for speed: Not set up for use as a base class
///  @li Represents 4x4 homogenous matrix as a 4x3 table with the last row
///      implicitly represented as [ 0, 0, 0, 1 ]

#ifndef INCLUDED_numeric_HomogeneousTransform_hh
#define INCLUDED_numeric_HomogeneousTransform_hh

#include <numeric/HomogeneousTransform.fwd.hh>
#include <numeric/constants.hh>
#include <numeric/IOTraits.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <iostream>

namespace numeric {

template  < class T >
class HomogeneousTransform {
public:
	typedef T value_type;

public:
	/// @brief Default constructor: axis aligned with the global coordinate frame
	/// and the point located at the origin.
	HomogeneousTransform() :
		xx_( 1.0 ),
		yx_( 0.0 ),
		zx_( 0.0 ),
		px_( 0.0 ),
		xy_( 0.0 ),
		yy_( 1.0 ),
		zy_( 0.0 ),
		py_( 0.0 ),
		xz_( 0.0 ),
		yz_( 0.0 ),
		zz_( 1.0 ),
		pz_( 0.0 )
	{
	}

	/// @brief Construct the coordinate frame from the points defined such that
	/// 1. the point is at p3,
	/// 2. the z axis points along p2 to p3.
	/// 3. the y axis is in the p1-p2-p3 plane
	/// 4. the x axis is the cross product of y and z.
	HomogeneousTransform(
		xyzVector< T > const & p1,
		xyzVector< T > const & p2,
		xyzVector< T > const & p3
	) :
		px_( p3.x() ),
		py_( p3.y() ),
		pz_( p3.z() )
	{
		// p1, p2 and p3 should never be the same.
		assert( p1.distance_squared( p2 ) > 1e-3 );
		assert( p1.distance_squared( p3 ) > 1e-3 );
		assert( p2.distance_squared( p3 ) > 1e-3 );


		xyzVector< T > v23 = p3 - p2;
		v23.normalize();
		zx_ = v23.x();
		zy_ = v23.y();
		zz_ = v23.z();

		xyzVector< T > v21 = p1 - p2;
		T dotprod = v21.dot( v23 );
		v21 -= dotprod * v23;
		v21.normalize();
		yx_ = v21.x();
		yy_ = v21.y();
		yz_ = v21.z();

		xyzVector< T > xaxis = v21.cross( v23 );
		xx_ = xaxis.x();
		xy_ = xaxis.y();
		xz_ = xaxis.z();

		assert( orthonormal() );

	}

	/// @brief Constructor from xyzVectors: trust that the axis are indeed orthoganol
	HomogeneousTransform(
		xyzVector< T > const & xaxis,
		xyzVector< T > const & yaxis,
		xyzVector< T > const & zaxis,
		xyzVector< T > const & point
	) :
		xx_( xaxis.x() ),
		yx_( yaxis.x() ),
		zx_( zaxis.x() ),
		px_( point.x() ),
		xy_( xaxis.y() ),
		yy_( yaxis.y() ),
		zy_( zaxis.y() ),
		py_( point.y() ),
		xz_( xaxis.z() ),
		yz_( yaxis.z() ),
		zz_( zaxis.z() ),
		pz_( point.z() )
	{
		assert( orthonormal() );
	}

	/// @brief Constructor from an xyzMatrix and xyzVectors: trust that the axis are indeed orthoganol
	/// @details Note that variable naming scheme (e.g. xy) is oposite of that in the xyzMatrix.
	HomogeneousTransform(
		xyzMatrix< T > const & axes,
		xyzVector< T > const & point
	) :
		xx_( axes.xx() ),
		yx_( axes.xy() ),
		zx_( axes.xz() ),
		px_( point.x() ),
		xy_( axes.yx() ),
		yy_( axes.yy() ),
		zy_( axes.yz() ),
		py_( point.y() ),
		xz_( axes.zx() ),
		yz_( axes.zy() ),
		zz_( axes.zz() ),
		pz_( point.z() )
	{
		assert( orthonormal() );
	}

	~HomogeneousTransform() {}

	HomogeneousTransform< T > const &
	operator = ( HomogeneousTransform< T > const & rhs ) {
		if ( this != & rhs ) {
			xx_ = rhs.xx_;
			xy_ = rhs.xy_;
			xz_ = rhs.xz_;
			yx_ = rhs.yx_;
			yy_ = rhs.yy_;
			yz_ = rhs.yz_;
			zx_ = rhs.zx_;
			zy_ = rhs.zy_;
			zz_ = rhs.zz_;
			px_ = rhs.px_;
			py_ = rhs.py_;
			pz_ = rhs.pz_;
		}
		return *this;
	}

public:
	//Accessors

	value_type xx() const { return xx_; }
	value_type xy() const { return xy_; }
	value_type xz() const { return xz_; }
	value_type yx() const { return yx_; }
	value_type yy() const { return yy_; }
	value_type yz() const { return yz_; }
	value_type zx() const { return zx_; }
	value_type zy() const { return zy_; }
	value_type zz() const { return zz_; }
	value_type px() const { return px_; }
	value_type py() const { return py_; }
	value_type pz() const { return pz_; }


	xyzVector< T >
	xaxis() const {
		return xyzVector< T >( xx_, xy_, xz_ );
	}

	xyzVector< T >
	yaxis() const {
		return xyzVector< T >( yx_, yy_, yz_ );
	}

	xyzVector< T >
	zaxis() const {
		return xyzVector< T >( zx_, zy_, zz_ );
	}

	xyzVector< T >
	point() const {
		return xyzVector< T >( px_, py_, pz_ );
	}

	xyzMatrix< T >
	rotation_matrix() const {

		// Entry naming is reversed from that of xyzMatrix,
		// see (xyzVector, xyzMatrix) constructor for specifc details.

		return xyzMatrix< T >::rows_constructor(
			xyzVector< T >( xx_, yx_, zx_ ),
			xyzVector< T >( xy_, yy_, zy_ ),
			xyzVector< T >( xz_, yz_, zz_ ));
	}


public:
	/// Mutators

	void
	set_identity_rotation() {
		xx_ = yy_ = zz_ = 1.0;
		xy_ = xz_ = yx_ = yz_ = zx_ = zy_ = 0.0;
	}

	void
	set_identity_transform() {
		px_ = py_ = pz_ = 0.0;
	}

	void
	set_identity() {
		set_identity_rotation();
		set_identity_transform();
	}

	void
	set_xaxis_rotation_deg( T angle ) {
		xyzMatrix< T > xrotmat = x_rotation_matrix_degrees( angle );
		HomogeneousTransform< T > temp( xrotmat, xyzVector< T >( 0.0 ) );
		*this = temp;
	}

	void
	set_yaxis_rotation_deg( T angle ) {
		xyzMatrix< T > yrotmat = y_rotation_matrix_degrees( angle );
		HomogeneousTransform< T > temp( yrotmat, xyzVector< T >( 0.0 ) );
		*this = temp;
	}

	void
	set_zaxis_rotation_deg( T angle ) {
		xyzMatrix< T > zrotmat = z_rotation_matrix_degrees( angle );
		HomogeneousTransform< T > temp( zrotmat, xyzVector< T >( 0.0 ) );
		*this = temp;
	}


	void
	set_xaxis_rotation_rad( T angle ) {
		xyzMatrix< T > xrotmat = x_rotation_matrix_radians( angle );
		HomogeneousTransform< T > temp( xrotmat, xyzVector< T >( 0.0 ) );
		*this = temp;
	}

	void
	set_yaxis_rotation_rad( T angle ) {
		xyzMatrix< T > yrotmat = y_rotation_matrix_radians( angle );
		HomogeneousTransform< T > temp( yrotmat, xyzVector< T >( 0.0 ) );
		*this = temp;
	}

	void
	set_zaxis_rotation_rad( T angle ) {
		xyzMatrix< T > zrotmat = z_rotation_matrix_radians( angle );
		HomogeneousTransform< T > temp( zrotmat, xyzVector< T >( 0.0 ) );
		*this = temp;
	}

	/// @brief Set this HT to describe a transformation along the three axes.
	/// This has the size effect of setting the axes to the global axes.
	void
	set_transform( xyzVector< T > const & t ) {
		set_identity_rotation();
		px_ = t.x();
		py_ = t.y();
		pz_ = t.z();
	}

	/// @brief Set the point that this HT is centered at.  This leaves
	/// the axes untouched.
	void
	set_point( xyzVector< T > const & p ) {
		px_ = p.x();
		py_ = p.y();
		pz_ = p.z();
	}

	/// @brief Right multiply this coordinate frame by the frame
	///   1  0  0  delta
	///   0  1  0    0
	///   0  0  1    0
	void
	walk_along_x( T delta ) {
		px_ += delta * xx_;
		py_ += delta * xy_;
		pz_ += delta * xz_;
	}

	/// @brief Right multiply this coordinate frame by the frame
	///   1  0  0    0
	///   0  1  0  delta
	///   0  0  1    0
	void
	walk_along_y( T delta ) {
		px_ += delta * yx_;
		py_ += delta * yy_;
		pz_ += delta * yz_;
	}

	/// @brief Right multiply this coordinate frame by the frame
	///   1  0  0    0
	///   0  1  0    0
	///   0  0  1  delta
	void
	walk_along_z( T delta ) {
		px_ += delta * zx_;
		py_ += delta * zy_;
		pz_ += delta * zz_;
	}

public:
	/// Multiplication

	/// @brief the main operation of a homogeneous transform: matrix multiplication.
	/// Note: matrix multiplication is transitive, so rotation/translation matrices
	/// may be pre-multiplied before being applied to a set of points.
	HomogeneousTransform< T >
	operator * ( HomogeneousTransform< T > const & rmat ) const {
		HomogeneousTransform< T > product;

		/// X axis
		product.xx_ = rmat.xx_ * xx_ + rmat.xy_ * yx_ + rmat.xz_ * zx_;
		product.xy_ = rmat.xx_ * xy_ + rmat.xy_ * yy_ + rmat.xz_ * zy_;
		product.xz_ = rmat.xx_ * xz_ + rmat.xy_ * yz_ + rmat.xz_ * zz_;

		/// Y axis
		product.yx_ = rmat.yx_ * xx_ + rmat.yy_ * yx_ + rmat.yz_ * zx_;
		product.yy_ = rmat.yx_ * xy_ + rmat.yy_ * yy_ + rmat.yz_ * zy_;
		product.yz_ = rmat.yx_ * xz_ + rmat.yy_ * yz_ + rmat.yz_ * zz_;

		/// Z axis
		product.zx_ = rmat.zx_ * xx_ + rmat.zy_ * yx_ + rmat.zz_ * zx_;
		product.zy_ = rmat.zx_ * xy_ + rmat.zy_ * yy_ + rmat.zz_ * zy_;
		product.zz_ = rmat.zx_ * xz_ + rmat.zy_ * yz_ + rmat.zz_ * zz_;

		/// Point
		product.px_ = rmat.px_ * xx_ + rmat.py_ * yx_ + rmat.pz_ * zx_ + px_;
		product.py_ = rmat.px_ * xy_ + rmat.py_ * yy_ + rmat.pz_ * zy_ + py_;
		product.pz_ = rmat.px_ * xz_ + rmat.py_ * yz_ + rmat.pz_ * zz_ + pz_;

		return product;
	}

	/// @brief Transform a point.  The input point is a location in this coordinate frame
	/// the output point is the location in global coordinate frame.
	xyzVector< T >
	operator * ( xyzVector< T > const & vect ) const {
		return xyzVector< T >(
			vect.x() * xx_ + vect.y() * yx_ + vect.z() * zx_ + px_,
			vect.x() * xy_ + vect.y() * yy_ + vect.z() * zy_ + py_,
			vect.x() * xz_ + vect.y() * yz_ + vect.z() * zz_ + pz_ );
	}

	/// @brief Invert this matrix.
	HomogeneousTransform< T >
	inverse() const {
		HomogeneousTransform< T > inv;
		inv.xx_ = xx_; inv.yx_ = xy_; inv.zx_ = xz_; inv.px_ = -( xx_ * px_ + xy_ * py_ + xz_ * pz_ );
		inv.xy_ = yx_; inv.yy_ = yy_; inv.zy_ = yz_; inv.py_ = -( yx_ * px_ + yy_ * py_ + yz_ * pz_ );
		inv.xz_ = zx_; inv.yz_ = zy_; inv.zz_ = zz_; inv.pz_ = -( zx_ * px_ + zy_ * py_ + zz_ * pz_ );
		return inv;
	}


	/// @brief Convert a point in the global coordinate system to a point in this coordinate frame.
	/// If this frame is F and the point is p, then this solves for x st:
	/// F x = p.  Equivalent to computing  (F.inverse() * p).point().
	xyzVector< T >
	to_local_coordinate( xyzVector< T > const & v ) const {
		xyzVector< T > local = v - point();
		return xyzVector< T >(
			xaxis().dot( local ),
			yaxis().dot( local ),
			zaxis().dot( local )
		);
	}

	/// @brief Return the three euler angles (in radians) that describe this HomogeneousTransform as the series
	/// of a Z axis rotation by the angle phi (returned in position 1 of the output vector), followed by
	/// an X axis rotation by the angle theta (returned in position 3 of the output vector), followed by another
	/// Z axis rotation by the angle psi (returned in position 2 of the output vector).
	/// This code is a modified version of Alex Z's code from r++.
	///
	/// @details
	/// The range of phi is [ -pi, pi ];
	/// The range of psi is [ -pi, pi ];
	/// The range of theta is [ 0, pi ];
	///
	/// The function pretends that this HomogeneousTransform is the result of these three transformations;
	/// if it were, then the rotation matrix would be
	///
	/// FIGURE 1:
	/// R = [
	///       cos(psi)cos(phi)-cos(theta)sin(phi)sin(psi)        cos(psi)sin(phi)+cos(theta)cos(phi)sin(psi)      sin(psi)sin(theta)
	///      -sin(psi)cos(phi)-cos(theta)sin(phi)cos(psi)       -sin(psi)sin(phi)+cos(theta)cos(phi)cos(psi)      cos(psi)sin(theta)
	///                   sin(theta)sin(phi)                                 -sin(theta)cos(phi)                        cos(theta)
	/// ]
	///
	/// where each axis above is represented as a ROW VECTOR (to be distinguished from the
	/// HomogeneousTransform's representation of axes as COLUMN VECTORS).
	///
	/// The zz_ coordinate gives away theta.
	/// Theta may be computed as acos( zz_ ), or, as Alex does it, asin( sqrt( 1 - zz^2))
	/// Since there is redundancy in theta, this function chooses a theta with a positive
	/// sin(theta): i.e. quadrants I and II.  Assuming we have a positive sin theta
	/// pushes phi and psi into conforming angles.
	///
	/// NOTE on theta: asin returns a value in the range [ -pi/2, pi/2 ], and we have artificially
	/// created a positive sin(theta), so we will get a asin( pos_sin_theta ), we have a value
	/// in the range [ 0, pi/2 ].  To convert this into the actual angle theta, we examine the zz sign.
	/// If zz is negative, we chose the quadrant II theta.
	/// That is, asin( pos_sin_theta) returned an angle, call it theta'.  Now, if cos( theta ) is negative,
	/// then we want to choose the positive x-axis rotation that's equivalent to -theta'.  To do so,
	/// we reflect q through the y axis (See figure 2 below) to get p and then measure theta as pi - theta'.
	///
	/// FIGURE 2:
	///
	///  II        |         I
	///            |
	///    p.      |      .q (cos(-theta'), std::abs(sin(theta')))
	///       .    |    .
	/// theta'( .  |  .  )  theta' = asin( std::abs(sin(theta))
	/// -----------------------
	///            |
	///            |
	///            |
	///  III       |        IV
	///            |
	///  The angle between the positive x axis and p is pi - theta'.
	///
	///
	///
	/// Since zx and zy contain only phi terms and a constant sin( theta ) term,
	/// phi is given by atan2( sin_phi, cos_phi ) = atan2( c*sin_phi, c*cos_phi ) = atan2( zx, -zy )
	/// for c positive and non-zero.  If sin_theta is zero, or very close to zero, we're at gimbal lock.
	///
	/// Moreover, since xz and yz contain only psi terms, psi may also be deduced using atan2.
	///
	/// There are 2 degenerate cases (gimbal lock)
	/// 1. theta close to 0  (North Pole singularity), or
	/// 2. theta close to pi (South Pole singularity)
	/// For these, we take: phi=acos(xx), theta = 0 (resp. Pi/2), psi = 0
	xyzVector< T >
	euler_angles_rad() const {
		xyzVector< T > euler;

		T const FLOAT_PRECISION( 1e-5 );

		// WARNING: Gimbal Lock!
		if ( zz_ >= 1 - FLOAT_PRECISION ) {
			euler(1) = -std::atan2( sin_cos_range( yx_), sin_cos_range( xx_ ) );
			euler(2) = 0.0;
			euler(3) = 0.0;
			return euler;
		}

		if ( zz_ <= -1 + FLOAT_PRECISION ) {
			euler(1) = -std::atan2( sin_cos_range( yx_), sin_cos_range( xx_ ) );
			euler(2) = 0.0;
			euler(3) = (T)numeric::constants::d::pi;
			return euler;
		}

		T pos_sin_theta = std::sqrt( 1 - zz_*zz_ ); // sin2theta = 1 - cos2theta.

		// two values are possible here: my convention is to use positive theta only.
		// corresponding theta between [0,pi/2] -> [0,90] since st > 0
		// and asin returns value between [-pi/2, pi/2]
		euler(3) = std::asin( pos_sin_theta );

		// decide whether the actual positive theta is between [pi/2, pi[ using the value of cos(theta)
		// which happens to be the matrix element zz_ (and is thus signed).
		if ( zz_ < 0 ) {
			euler(3) = numeric::constants::d::pi - euler(3);
		}

		//euler(1) = std::atan2(  -UU(1,3), UU(2,3) );   // between -Pi and Pi -> [-180,180]
		//euler(2) = std::atan2( UU(3,1), UU(3,2) );     // between -Pi and Pi -> [-180, 180]

		// this is atan( sin_phi * c, cos_phi * c  ) as opposed to Alex's atan( -sin_phi * c, -cos_phi * c ).
		euler(1) = std::atan2( zx_, -zy_ );
		euler(2) = std::atan2( xz_, yz_ );


		return euler;
	}

	xyzVector< T >
	euler_angles_deg() const {
		return numeric::constants::d::radians_to_degrees * euler_angles_rad();
	}

	/// @brief Construct the coordinate frame from three euler angles that describe the frame.
	/// Keep the point fixed.  See the description for euler_angles_rad() to understand
	/// the Z-X-Z transformation convention.
	void
	from_euler_angles_rad( xyzVector< T > const & euler ) {
		/*
		HomogeneousTransform< T > zrot1, xrot2, zrot3;
		zrot1.set_zaxis_rotation_rad( euler( 1 ) );
		xrot2.set_xaxis_rotation_rad( euler( 3 ) );
		zrot3.set_zaxis_rotation_rad( euler( 2 ) );
		T const px( px_), py( py_ ), pz( pz_ );
		(*this) = zrot1 * xrot2 * zrot3;
		px_ = px; py_ = py; pz_ = pz;
		*/

		T const phi( euler( 1 ) ), psi( euler( 2 ) ), theta( euler( 3 ) );

		T const cos_phi( std::cos( phi ) ),    sin_phi( std::sin( phi ) );
		T const cos_psi( std::cos( psi ) ),    sin_psi( std::sin( psi ) );
		T const cos_theta( std::cos( theta )), sin_theta( std::sin( theta ) );


		xx_ =  cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;  xy_ =  cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;  xz_ =  sin_psi * sin_theta;
		yx_ = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi;  yy_ = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi;  yz_ =  cos_psi * sin_theta;
		zx_ =                sin_theta * sin_phi;                  zy_ =                     -sin_theta * cos_phi;            zz_ =        cos_theta;

	}

	void
	from_euler_angles_deg( xyzVector< T > const & euler ) {
		xyzVector< T > euler_rad( euler );
		euler_rad *= numeric::constants::d::degrees_to_radians;
		from_euler_angles_rad( euler_rad );
	}

public:
	std::ostream &
	show_stream( std::ostream & stream = std::cout) const
	{
		// Types
		using std::setw;
		typedef IOTraits< T >  Traits;

		// Save current stream state and set persistent state
		std::ios_base::fmtflags const old_flags = stream.flags();
		int const old_precision = stream.precision( Traits::precision() );
		stream << std::right << std::showpoint << std::uppercase;

		// Output xyzVector
		int const w = Traits::width();
		stream << setw( w ) << xx() << ' ' << setw( w ) << yx() << ' ' << setw( w ) << zx() << setw( w ) << px() << std::endl;
		stream << setw( w ) << xy() << ' ' << setw( w ) << yy() << ' ' << setw( w ) << zy() << setw( w ) << py() << std::endl;
		stream << setw( w ) << xz() << ' ' << setw( w ) << yz() << ' ' << setw( w ) << zz() << setw( w ) << pz() << std::endl;


		// Restore previous stream state
		stream.precision( old_precision );
		stream.flags( old_flags );

		return stream;
	}

	void
	show(std::ostream & stream = std::cout) const
	{
		stream << show_stream( stream );
	}

private:

	bool
	orthonormal() const {
		return orthoganol() && normal();
	}

	/// @brief Are the axis orthoganol
	bool
	orthoganol() const {
		if ( xx_*yx_ + xy_*yy_ + xz_*yz_ > 1e-6 ) return false;
		if ( xx_*zx_ + xy_*zy_ + xz_*zz_ > 1e-6 ) return false;
		if ( yx_*zx_ + yy_*zy_ + yz_*zz_ > 1e-6 ) return false;
		return true;
	}

	/// @brief Are the axis of unit length?
	bool
	normal() const {
		if ( std::abs( xx_*xx_ + xy_*xy_ + xz_*xz_ -1 ) > 1e-6 ) return false;
		if ( std::abs( yx_*yx_ + yy_*yy_ + yz_*yz_ -1 ) > 1e-6 ) return false;
		if ( std::abs( zx_*zx_ + zy_*zy_ + zz_*zz_ -1 ) > 1e-6 ) return false;
		return true;
	}

private:
	/// Variables below are allocated in an intentionally visually pleasing order
	/// for a code reader, not necessarily for performance.
	/// The homogenous matrix is a 4 x 3 matrix, with a pseudo row-major ordering
	/// of data.

	/// Naming scheme:
	/// Column 1 is the x axis as a column vector
	/// Column 2 is the y axis as a column vector
	/// Column 3 is the z axis as a column vector
	/// Column 4 is the location of the point, as a column vector.

	/// Row 1 is the x coordinate of the three axes and the point
	/// Row 2 is the y coordinate of the three axes and the point
	/// Row 3 is the z coordinate of the three axes and the point
	/// Row 4 is not explicitly represented, but is always [ 0, 0, 0, 1 ]

	/// For axis Q, the R component is named qr_; e.g. The Z component of the Y axis is yz_;

	/// Note: This naming scheme is different from the one that Stuart uses in the
	/// xyzMatrix class.

	value_type xx_, yx_, zx_, px_;
	value_type xy_, yy_, zy_, py_;
	value_type xz_, yz_, zz_, pz_;


};

template< typename T >
std::ostream &
operator <<( std::ostream & stream, HomogeneousTransform< T > const & ht )
{
	return ht.show_stream( stream );
}

// PyRosetta WorkAround
class HomogeneousTransform_Double : public HomogeneousTransform<double>
{
public:
	HomogeneousTransform_Double() {}

	HomogeneousTransform_Double(
		xyzVector< double > const & p1,
		xyzVector< double > const & p2,
		xyzVector< double > const & p3
	) : HomogeneousTransform<double>(p1, p2, p3) {}

	HomogeneousTransform_Double(
		xyzVector< double > const & xaxis,
		xyzVector< double > const & yaxis,
		xyzVector< double > const & zaxis,
		xyzVector< double > const & point) : HomogeneousTransform<double>(xaxis,yaxis,zaxis,point) {}

	HomogeneousTransform_Double(
		xyzMatrix< double > const & axes,
		xyzVector< double > const & point) : HomogeneousTransform<double>(axes, point) {}
};

inline std::ostream &
operator << ( std::ostream & stream, HomogeneousTransform< double > const & ht )
{
	return ht.show_stream( stream );
}

}

#endif
