// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/EulerAngles.hh
/// @brief  EulerAngles container derived from xyzVector.hh
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_EulerAngles_hh
#define INCLUDED_numeric_EulerAngles_hh

// Unit headers
#include <numeric/EulerAngles.fwd.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/types.hh>

// C++ headers
#include <cmath>

namespace numeric {

/// @brief Euler angles 3-D orientation representation
///
/// @remarks
/// The three euler angles (in radians) that describing a rotation operation
/// of a Z axis rotation by the angle phi (position 1), followed by
/// an X axis rotation by the angle theta (position 3), followed by another
/// Z axis rotation by the angle psi (position 2).
/// this->code is a modified version of Alex Z's code from r++.
///
/// @details
/// The range of phi is ( -pi, pi ];
/// The range of psi is ( -pi, pi ];
/// The range of theta is [ 0, pi ];
template< typename T >
class EulerAngles : public xyzVector<T>
{
	typedef xyzVector<T> Vector;
	typedef xyzMatrix<T> Matrix;
	typedef T Value;

public:

	/// @brief Default constructor
	/// @note  Values are initialized to zero
	inline
	EulerAngles() : xyzVector<T>(0.)
	{}

	/// @brief Copy constructor
	inline
	EulerAngles( xyzVector< T > const & v ) : xyzVector<T>(v)
	{}

	/// @brief Triple value constructor
	inline
	EulerAngles(
		Value const & phi,
		Value const & psi,
		Value const & theta
	) : xyzVector<T>( phi, psi, theta )
	{}

	/// @brief Construct from rotation matrix.
	inline EulerAngles(xyzMatrix<T> rotation_matrix) : xyzVector<T>()
	{
		from_rotation_matrix(rotation_matrix);
	}

	/// @brief Destructor
	inline
	~EulerAngles()
	{}


	/// The equivalent rotation matrix representation of the euler angles would be:
	///
	/// FIGURE 1:
	/// R = [
	///       cos(psi)cos(phi)-cos(theta)sin(phi)sin(psi)        cos(psi)sin(phi)+cos(theta)cos(phi)sin(psi)      sin(psi)sin(theta)
	///      -sin(psi)cos(phi)-cos(theta)sin(phi)cos(psi)       -sin(psi)sin(phi)+cos(theta)cos(phi)cos(psi)      cos(psi)sin(theta)
	///                   sin(theta)sin(phi)                                 -sin(theta)cos(phi)                        cos(theta)
	/// ]
	///
	/// The zz_ coordinate gives away theta.
	/// Theta may be computed as acos( zz_ ), or, as Alex does it, asin( sqrt( 1 - zz^2))
	/// Since there is redundancy in theta, this->function chooses a theta with a positive
	/// sin(theta): i.e. quadrants I and II.  Assuming we have a positive sin theta
	/// pushes phi and psi into conforming angles.
	///
	/// NOTE on theta: asin returns a value in the range [ -pi/2, pi/2 ], and we have artificially
	/// created a positive sin(theta), so we will get a asin( pos_sin_theta ), we have a value
	/// in the range [ 0, pi/2 ].  To convert this->into the actual angle theta, we examine the zz sign.
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

	void
	from_rotation_matrix(xyzMatrix<T> matrix)
	{
		T const FLOAT_PRECISION( 1e-5 );

		// WARNING: Gimbal Lock!
		if ( matrix.zz() >= 1 - FLOAT_PRECISION ) {
			phi(std::acos( sin_cos_range( matrix.xx() ) ));
			psi(0.0);
			theta(0.0);
			return;
		}

		if ( matrix.zz() <= -1 + FLOAT_PRECISION ) {
			phi(std::acos( sin_cos_range( matrix.xx() )));
			psi(0.0);
			theta(numeric::constants::d::pi);
			return;
		}

		T pos_sin_theta = std::sqrt( 1 - matrix.zz() * matrix.zz()); // sin2theta = 1 - cos2theta.

		// two values are possible here: my convention is to use positive theta only.
		// corresponding theta between [0,pi/2] -> [0,90] since st > 0
		// and asin returns value between [-pi/2, pi/2]
		theta(std::asin( pos_sin_theta ));

		// decide whether the actual positive theta is between [pi/2, pi[ using the value of cos(theta)
		// which happens to be the matrix element zz_ (and is thus signed).
		if ( matrix.zz() < 0 ) {
			theta(numeric::constants::d::pi - theta());
		}

		//euler(1) = std::atan2(  -UU(1,3), UU(2,3) );   // between -Pi and Pi -> [-180,180]
		//euler(2) = std::atan2( UU(3,1), UU(3,2) );     // between -Pi and Pi -> [-180, 180]

		// this->is atan( sin_phi * c, cos_phi * c  ) as opposed to Alex's atan( -sin_phi * c, -cos_phi * c ).
		phi(std::atan2( matrix.xz(), -matrix.yz()));
		psi(std::atan2( matrix.zx(), matrix.zy()));

		// Check and correct -pi to pi
		if ( phi() + numeric::constants::d::pi <= FLOAT_PRECISION ) {
			phi(numeric::constants::d::pi);
		}
		if ( psi() + numeric::constants::d::pi <= FLOAT_PRECISION ) {
			psi(numeric::constants::d::pi);
		}
	}

	/// @brief Construct rotation matrix from three euler angles that describe the frame.
	/// See the description for from_rotation_matrix to understand
	/// the Z-X-Z transformation convention.
	xyzMatrix<T>
	to_rotation_matrix()
	{
		xyzMatrix<T> m;

		T const cos_phi( std::cos( phi() ) ),    sin_phi( std::sin( phi() ) );
		T const cos_psi( std::cos( psi() ) ),    sin_psi( std::sin( psi() ) );
		T const cos_theta( std::cos( theta() )), sin_theta( std::sin( theta() ) );

		m.xx(cos_psi * cos_phi - cos_theta * sin_phi * sin_psi );
		m.yx( cos_psi * sin_phi + cos_theta * cos_phi * sin_psi );
		m.zx( sin_psi * sin_theta );

		m.xy( -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi );
		m.yy( -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi );
		m.zy(  cos_psi * sin_theta );

		m.xz( sin_theta * sin_phi );
		m.yz( -sin_theta * cos_phi );
		m.zz( cos_theta );

		return m;
	}

	/// @brief Value phi in radians
	inline
	Value &
	phi()
	{
		return Vector::x();
	}

	/// @brief Set value phi in radians
	inline
	void
	phi(Value const & value )
	{
		Vector::x(value);
	}

	/// @brief Value phi in radians
	inline
	Value &
	phi_radians()
	{
		return Vector::x();
	}

	/// @brief Set value phi in radians
	inline
	void
	phi_radians(Value const & value )
	{
		Vector::x(value);
	}

	/// @brief Value phi in degrees
	inline
	Value
	phi_degrees()
	{
		return numeric::constants::d::radians_to_degrees * Vector::x();
	}

	/// @brief Set value phi in degrees
	inline
	void
	phi_degrees(Value const & value )
	{
		Vector::x(value * numeric::constants::d::degrees_to_radians);
	}

	/// @brief Value psi in radians
	inline
	Value &
	psi()
	{
		return Vector::y();
	}

	/// @brief Set value psi in radians
	inline
	void
	psi(Value const & value )
	{
		Vector::y(value);
	}

	/// @brief Value psi in radians
	inline
	Value &
	psi_radians()
	{
		return Vector::y();
	}

	/// @brief Set value psi in radians
	inline
	void
	psi_radians(Value const & value )
	{
		Vector::y(value);
	}

	/// @brief Value psi in degrees
	inline
	Value
	psi_degrees()
	{
		return numeric::constants::d::radians_to_degrees * Vector::y();
	}

	/// @brief Set value psi in degrees
	inline
	void
	psi_degrees(Value const & value )
	{
		Vector::y(value * numeric::constants::d::degrees_to_radians);
	}

	/// @brief Value theta in radians
	inline
	Value &
	theta()
	{
		return Vector::z();
	}

	/// @brief Set value theta in radians
	inline
	void
	theta(Value const & value )
	{
		Vector::z(value);
	}

	/// @brief Value theta in radians
	inline
	Value &
	theta_radians()
	{
		return Vector::z();
	}

	/// @brief Set value theta in radians
	inline
	void
	theta_radians(Value const & value )
	{
		Vector::z(value);
	}

	/// @brief Value theta in degrees
	inline
	Value
	theta_degrees()
	{
		return numeric::constants::d::radians_to_degrees * Vector::z();
	}

	/// @brief Set value theta in degrees
	inline
	void
	theta_degrees(Value const & value )
	{
		Vector::z(value * numeric::constants::d::degrees_to_radians);
	}

	/// @brief Get angular distance between two sets of Euler Angles.
	static
	T
	angular_distance_between(EulerAngles<T> a1, EulerAngles<T> a2)
	{
		xyzMatrix<T> rotation = a1.to_rotation_matrix().transposed() * a2.to_rotation_matrix();

		return rotation_angle(rotation);
	}

	// @brief Static constructor from degrees
	static
	EulerAngles<T>
	from_degrees(T phi, T psi, T theta)
	{
		return EulerAngles<T>(
			conversions::radians(phi),
			conversions::radians(psi),
			conversions::radians(theta)
		);
	}

	// @brief Static constructor from degrees
	static
	EulerAngles<T>
	from_degrees(xyzVector<T> vector)
	{
		return EulerAngles<T>(
			conversions::radians(vector.x()),
			conversions::radians(vector.y()),
			conversions::radians(vector.z())
		);
	}

	// @brief Static constructor from radians
	static
	EulerAngles<T>
	from_radians(T phi, T psi, T theta)
	{
		return EulerAngles<T>(phi, psi, theta);
	}

	// @brief Static constructor from degrees
	static
	EulerAngles<T>
	from_radians(xyzVector<T> vector)
	{
		return EulerAngles<T>(vector);
	}
};

} // namespace numeric

#endif // INCLUDED_numeric_EulerAngles_hh
