// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/BodyPosition.hh
/// @brief  Rigid body 3-D position/transform
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_BodyPosition_hh
#define INCLUDED_numeric_BodyPosition_hh


// Unit headers
#include <numeric/BodyPosition.fwd.hh>

// Package headers
#include <numeric/Quaternion.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>

// C++ headers
#include <utility/assert.hh>
#include <cmath>


namespace numeric {


/// @brief Rigid body 3-D position/transform
///
/// @remarks
///  @li Dual personality: rotation matrix and quaternion representations
///  @li Rotation matrix is faster for transforming points
///  @li Quaternion is faster for successive frame transforms and easier to keep normalized
///  @li For quaternion q = ( w, x, y, z ) the equivalent rotation matrix is:
/// 	     1-2( y^2 + z^2 )  2( x y - w z )    2( x z + w y )
/// 	     2( x y + w z )    1-2( x^2 + z^2 )  2( y z - w x )
/// 	     2( x z - w y )    2( y z + w x )    1-2( x^2 + y^2 )
///  @li Quaternions q and -q have the same rotation matrix
///  @li If an odd multiple of 2*pi is added to an angle the quaternion q is changed to -q
///  @li Invariant: The rotation matrix or quaternion or both must be "fresh" before and
///      after every member function call
template< typename T >
class BodyPosition
{


public: // Types


	// Project style
	typedef  T          Value;
	typedef  T &        Reference;
	typedef  T const &  ConstReference;
	typedef  T *        Pointer;
	typedef  T const *  ConstPointer;

	// STL/boost style
	typedef  T          value_type;
	typedef  T &        reference;
	typedef  T const &  const_reference;
	typedef  T *        pointer;
	typedef  T const *  const_pointer;

	// Traits
	typedef  NumericTraits< T >  Traits;

	// Math
	typedef  xyzVector< T >  Axis;
	typedef  xyzVector< T >  Point;
	typedef  xyzVector< T >  Vector;
	typedef  xyzVector< T >  Translation;
	typedef  xyzMatrix< T >  Rotation;
	typedef  numeric::Quaternion< T >  Quaternion;


public: // Creation


	/// @brief Default constructor
	inline
	BodyPosition() :
		R_( Rotation::identity() ),
		R_fresh_( true ),
		q_fresh_( true ),
		t_( T( 0 ) )
	{}


	/// @brief Rotation + Translation constructor
	inline
	BodyPosition(
		Rotation const & R_a,
		Translation const & t_a
	) :
		R_( R_a ),
		R_fresh_( true ),
		q_fresh_( false ),
		t_( t_a )
	{}


	/// @brief Rotation constructor
	inline
	BodyPosition( Rotation const & R_a ) :
		R_( R_a ),
		R_fresh_( true ),
		q_fresh_( false ),
		t_( T( 0 ) )
	{}


	/// @brief Quaternion + Translation constructor
	inline
	BodyPosition(
		Quaternion const & q_a,
		Translation const & t_a
	) :
		R_fresh_( false ),
		q_( q_a ),
		q_fresh_( true ),
		t_( t_a )
	{}


	/// @brief Quaternion constructor
	inline
	BodyPosition( Quaternion const & q_a ) :
		R_fresh_( false ),
		q_( q_a ),
		q_fresh_( true ),
		t_( T( 0 ) )
	{}


	/// @brief Point + Vector constructor
	/// @note  Constructs the transform that takes a point p and unit vectors x and u
	///        bound to p and separated by an angle on (0,pi) from position 1 to 2
	inline
	BodyPosition(
		Point const &,
		Vector const & x1,
		Vector const & u1,
		Point const &,
		Vector const & x2,
		Vector const & u2,
		Value const & cos,
		Value const & sin,
		Value const & csc
	) :
		R_fresh_( true ),
		q_fresh_( false ),
		t_( T( 0 ) )
	{
		assert( x1.is_normalized( Traits::length_tolerance() ) );
		assert( u1.is_normalized( Traits::length_tolerance() ) );
		assert( x2.is_normalized( Traits::length_tolerance() ) );
		assert( u2.is_normalized( Traits::length_tolerance() ) );
		assert( ( cos > T( -1 ) ) && ( cos < T( 1 ) ) );
		assert( ( sin > T( 0 ) ) && ( sin <= T( 1 ) ) );
//		assert( equal( square( cos ) + square( sin ), T( 1 ), square( Traits::sin_cos_tolerance() ) ) );
//		assert( equal( dot( x1, u1 ), cos, Traits::sin_cos_tolerance() ) );
//		assert( equal( dot( x2, u2 ), cos, Traits::sin_cos_tolerance() ) );
		Vector y1( x1 ); y1 *= -cos += u1 *= csc;
//		R_ = Rotation::rows(
//			x1,
//			y1,
//			cross( x1, y1 )
//		);
	}


	/// @brief Identity named constructor
	inline
	static
	BodyPosition
	identity()
	{
		return BodyPosition(); // Default constructor gives identity
	}


	/// @brief Copy constructor
	inline
	BodyPosition( BodyPosition const & bp ) :
		R_( bp.R_ ),
		R_fresh_( bp.R_fresh_ ),
		q_( bp.q_ ),
		q_fresh_( bp.q_fresh_ ),
		t_( bp.t_ )
	{}


	/// @brief Destructor
	inline
	~BodyPosition()
	{}


public: // assignment

	/// @brief Copy assignment
	inline
	BodyPosition &
	operator =( BodyPosition const & bp )
	{
		if ( this != &bp ) {
			R_ = bp.R_;
			R_fresh_ = bp.R_fresh_;
			q_ = bp.q_;
			q_fresh_ = bp.q_fresh_;
			t_ = bp.t_;
		}
		return *this;
	}


public: // Properties


	/// @brief Rotation
	inline
	Rotation const &
	R() const
	{
		if ( ! R_fresh_ ) R_refresh();
		return R_;
	}


	/// @brief Quaternion
	inline
	Quaternion const &
	q() const
	{
		if ( ! q_fresh_ ) q_refresh();
		return q_;
	}


	/// @brief Translation
	inline
	Translation const &
	t() const
	{
		return t_;
	}


	/// @brief Principal angle of rotation (on [0,2*pi])
	inline
	Value
	angle() const
	{
		if ( ! q_fresh_ ) q_refresh();
		return q_.angle();
	}


	/// @brief Axis of Rotation unit vector (direction for angle on [0,2*pi])
	inline
	Axis
	axis() const
	{
		if ( ! q_fresh_ ) q_refresh();
		return q_.axis(); // Returns zero vector if angle is zero
	}


	/// @brief Axis of rotation unit vector: Passed vector (direction for angle on [0,2*pi])
	inline
	Axis &
	axis( Axis & u ) const
	{
		if ( ! q_fresh_ ) q_refresh();
		return q_.axis( u ); // Returns zero vector if angle is zero
	}


public: // Methods: transforms


	/// @brief Transform a Point
	inline
	BodyPosition const &
	operator ()( Point & p ) const
	{
		if ( ! R_fresh_ ) R_refresh(); // Rotation matrix is more efficient than quaternion for this
		inplace_product( R_, p ) += t_;
		return *this;
	}


	/// @brief Transform a Point
	inline
	BodyPosition const &
	transform( Point & p ) const
	{
		if ( ! R_fresh_ ) R_refresh(); // Rotation matrix is more efficient than quaternion for this
		inplace_product( R_, p ) += t_;
		return *this;
	}


	/// @brief Inverse transform a Point
	inline
	BodyPosition const &
	inverse_transform( Point & p ) const
	{
		if ( ! R_fresh_ ) R_refresh(); // Rotation matrix is more efficient than quaternion for this
		inplace_transpose_product( R_, p -= t_ );
		return *this;
	}


public: // Methods: modifiers


	/// @brief Normalize
	inline
	BodyPosition &
	normalize()
	{
		if ( ! q_fresh_ ) q_refresh();
		q_.normalize();
		R_fresh_ = false;
		return *this;
	}


	/// @brief Normalize if magnitude squared error exceeds tolerance
	inline
	BodyPosition &
	normalize_if_needed( Value const & tol = Traits::quaternion_tolerance() )
	{
		if ( ! q_fresh_ ) q_refresh();
		q_.normalize_if_needed( tol );
		R_fresh_ = false;
		return *this;
	}


	/// @brief Identity
	inline
	BodyPosition &
	to_identity()
	{
		R_.to_identity();
		R_fresh_ = true;
		q_.to_identity();
		q_fresh_ = true;
		t_ = T( 0 );
		return *this;
	}


	/// @brief Invert
	inline
	BodyPosition &
	invert()
	{
		if ( t_ != T( 0 ) ) { // t -> -( R^T * t )
			if ( ! R_fresh_ ) R_refresh();
			R_.transpose(); // R^-1 == R^T since R is orthogonal
			inplace_product( R_, t_ );
			t_ *= -1;
		} else {
			if ( R_fresh_ ) R_.transpose(); // R^-1 == R^T since R is orthogonal
		}
		if ( q_fresh_ ) q_.invert();
		return *this;
	}


	/// @brief Apply a successive BodyPosition transformation
	inline
	BodyPosition &
	apply( BodyPosition const & p, bool const precise = true )
	{
		return left_transform_by( p, precise ); // OR SHOULD THIS BE RIGHT?
	}


	/// @brief Left transform by a BodyPosition
	inline
	BodyPosition &
	left_transform_by( BodyPosition const & p, bool const precise = true )
	{
		if ( ! q_fresh_ ) q_refresh();
		if ( ! p.q_fresh_ ) p.q_refresh();
		q_.left_multiply_by( p.q_, precise );
		R_fresh_ = false;
		p.transform( t_ );
		return *this;
	}


	/// @brief Right transform by a BodyPosition
	inline
	BodyPosition &
	right_transform_by( BodyPosition const & p, bool const precise = true )
	{
		if ( ! q_fresh_ ) q_refresh();
		if ( ! p.q_fresh_ ) p.q_refresh();
		q_.right_multiply_by( p.q_, precise );
		R_fresh_ = false;
		t_ += transformed( p.t_ );
		return *this;
	}


	/// @brief Left transform by the inverse of a BodyPosition
	inline
	BodyPosition &
	left_transform_by_inverse_of( BodyPosition const & p, bool const precise = true )
	{
		if ( ! q_fresh_ ) q_refresh();
		if ( ! p.q_fresh_ ) p.q_refresh();
		q_.left_multiply_by_inverse_of( p.q_, precise );
		R_fresh_ = false;
		p.inverse_transform( t_ );
		return *this;
	}


	/// @brief Right transform by the inverse of a BodyPosition
	inline
	BodyPosition &
	right_transform_by_inverse_of( BodyPosition const & p, bool const precise = true )
	{
		if ( ! q_fresh_ ) q_refresh();
		if ( ! p.q_fresh_ ) p.q_refresh();
		q_.right_multiply_by_inverse_of( p.q_, precise );
		R_fresh_ = false;
		t_ = transform( p.inverse_translation() );
		return *this;
	}


public: // Methods: generators


	/// @brief Inverse
	inline
	BodyPosition
	inverse() const
	{
		return BodyPosition( *this ).invert();
	}


	/// @brief Transformed Point
	inline
	Point
	transformed( Point const & p ) const
	{
		if ( ! R_fresh_ ) R_refresh(); // Rotation matrix is more efficient than quaternion for this
		return product( R_, p ) += t_;
	}


	/// @brief Inverse transformed Point
	inline
	Point
	inverse_transformed( Point const & p ) const
	{
		if ( ! R_fresh_ ) R_refresh(); // Rotation matrix is more efficient than quaternion for this
		return transpose_product( R_, p - t_ );
	}


	/// @brief Inverse translation
	inline
	Translation
	inverse_translation() const
	{
		if ( ! R_fresh_ ) R_refresh(); // Rotation matrix is more efficient than quaternion for this
		return transpose_product( R_, t_ ) *= -1; // -R^T * t
	}


//	/// @brief BodyPosition * BodyPosition
//	friend
//	inline
//	BodyPosition
//	operator *( BodyPosition const & p2, BodyPosition const & p1 )
//	{
//		return BodyPosition(
//		).normalize();
//	}
//
//
//	/// @brief Product: BodyPosition * BodyPosition
//	friend
//	inline
//	BodyPosition
//	product( BodyPosition const & p2, BodyPosition const & p1, bool const precise = true )
//	{
//		if ( precise ) { // Return normalized product
//			return BodyPosition(
//			).normalize();
//		} else { // Return product
//			return BodyPosition(
//			);
//		}
//	}


	/// @brief Identity BodyPosition for expressions
	/// @note  Default and identity() named constructors can be faster for construction
	/// @note  Can be safely used in construction of global objects
	inline
	static
	BodyPosition const &
	I()
	{
		static BodyPosition const I_; // Default constructor gives identity
		return I_;
	}


public: // Comparison


	/// @brief BodyPosition == BodyPosition
	friend
	inline
	bool
	operator ==( BodyPosition const & p1, BodyPosition const & p2 )
	{
		if ( p1.q_fresh_ && p2.q_fresh_ ) { // Use quaternions
			return ( ( p1.t_ == p2.t_ ) && ( p1.q_ == p2.q_ ) );
		} else if ( p1.R_fresh_ && p2.R_fresh_ ) { // Use rotations
			return ( ( p1.t_ == p2.t_ ) && ( p1.R_ == p2.R_ ) );
		} else { // Use quaternions
			if ( ! p1.q_fresh_ ) p1.q_refresh();
			if ( ! p2.q_fresh_ ) p2.q_refresh();
			return ( ( p1.t_ == p2.t_ ) && ( p1.q_ == p2.q_ ) );
		}
	}


	/// @brief BodyPosition != BodyPosition
	friend
	inline
	bool
	operator !=( BodyPosition const & p1, BodyPosition const & p2 )
	{
		return !( p1 == p2 );
	}


private: // Functions


	/// @brief Refresh rotation matrix from Quaternion
	inline
	BodyPosition const &
	R_refresh( bool const precise = true ) const
	{
		assert( q_fresh_ );

		if ( precise ) q_.normalize();

		Value const q_xs( q_.x_ * q_.x_ );
		Value const q_ys( q_.y_ * q_.y_ );
		Value const q_zs( q_.z_ * q_.z_ );
		Value const q_wx( q_.w_ * q_.x_ );
		Value const q_wy( q_.w_ * q_.y_ );
		Value const q_wz( q_.w_ * q_.z_ );
		Value const q_xy( q_.x_ * q_.y_ );
		Value const q_xz( q_.x_ * q_.z_ );
		Value const q_yz( q_.y_ * q_.z_ );

		R_.xx() = T( 1 ) - ( T( 2 ) * ( q_ys + q_zs ) );
		R_.xy() = T( 2 ) * ( q_xy - q_wz );
		R_.yx() = T( 2 ) * ( q_xy + q_wz );
		R_.yy() = T( 1 ) - ( T( 2 ) * ( q_xs + q_zs ) );
		R_.yz() = T( 2 ) * ( q_yz - q_wx );
		R_.zy() = T( 2 ) * ( q_yz + q_wx );
		R_.zz() = T( 1 ) - ( T( 2 ) * ( q_xs + q_ys ) );
		R_.zx() = T( 2 ) * ( q_xz - q_wy );
		R_.xz() = T( 2 ) * ( q_xz + q_wy );

		R_fresh_ = true;

		return *this;
	}


	/// @brief Refresh Quaternion from rotation matrix
	inline
	BodyPosition const &
	q_refresh( bool const precise = true ) const
	{
		assert( R_fresh_ );

		Value const d( T( 1 ) + R_.xx() + R_.yy() + R_.zz() ); // 1 + trace == 4(w^2)
		if ( d > T( 0 ) ) {
			Value const two_w( std::sqrt( d ) ); // 2w
			Value const fac( T( 0.5 ) / two_w ); // 1/4W
			q_.w_ = T( 0.5 ) * two_w;
			q_.x_ = ( R_.zy() - R_.yz() ) * fac;
			q_.y_ = ( R_.xz() - R_.zx() ) * fac;
			q_.z_ = ( R_.yx() - R_.xy() ) * fac;
		} else if ( ( R_.xx() >= R_.yy() ) && ( R_.xx() >= R_.zz() ) ) {
			Value const s( T( 2 ) * std::sqrt( T( 1 ) + R_.xx() - R_.yy() - R_.zz() ) ); // 4x
			Value const s_inv( T( 1 ) / s ); // 1/4x
			q_.w_ = ( R_.zy() - R_.yz() ) * s_inv;
			q_.x_ = T( 0.25 ) * s;
			q_.y_ = ( R_.xy() + R_.yx() ) * s_inv;
			q_.z_ = ( R_.xz() + R_.zx() ) * s_inv;
		} else if ( R_.yy() >= R_.zz() ) {
			Value const s( T( 2 ) * std::sqrt( T( 1 ) + R_.yy() - R_.xx() - R_.zz() ) ); // 4y
			Value const s_inv( T( 1 ) / s ); // 1/4y
			q_.w_ = ( R_.xz() - R_.zx() ) * s_inv;
			q_.x_ = ( R_.xy() + R_.yx() ) * s_inv;
			q_.y_ = T( 0.25 ) * s;
			q_.z_ = ( R_.yz() + R_.zy() ) * s_inv;
		} else {
			Value const s( T( 2 ) * std::sqrt( T( 1 ) + R_.zz() - R_.xx() - R_.yy() ) ); // 4z
			Value const s_inv( T( 1 ) / s ); // 1/4z
			q_.w_ = ( R_.yx() - R_.xy() ) * s_inv;
			q_.x_ = ( R_.xz() + R_.zx() ) * s_inv;
			q_.y_ = ( R_.yz() + R_.zy() ) * s_inv;
			q_.z_ = T( 0.25 ) * s;
		}

		q_fresh_ = true;

		if ( precise ) {
			q_.normalize();
			R_refresh();
		}

		return *this;
	}


private: // Fields


	/// @brief Rotation matrix
	mutable Rotation R_;

	/// @brief Rotation matrix status
	mutable bool R_fresh_;

	/// @brief Quaternion
	mutable Quaternion q_;

	/// @brief Quaternion status
	mutable bool q_fresh_;

	/// @brief Translation
	Translation t_;


}; // BodyPosition


/// @brief BodyPosition == BodyPosition
template< typename T >
bool
operator ==( BodyPosition< T > const & p1, BodyPosition< T > const & p2 );


/// @brief BodyPosition != BodyPosition
template< typename T >
bool
operator !=( BodyPosition< T > const & p1, BodyPosition< T > const & p2 );


} // namespace numeric


#endif // INCLUDED_numeric_BodyPosition_HH
