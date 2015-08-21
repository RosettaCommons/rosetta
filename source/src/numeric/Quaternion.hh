// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/Quaternion.hh
/// @brief  Unit quaternion 3-D orientation representation
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_Quaternion_hh
#define INCLUDED_numeric_Quaternion_hh


// Unit headers
#include <numeric/Quaternion.fwd.hh>

// Package headers
#include <numeric/NumericTraits.hh>
#include <numeric/xyzVector.hh>
#include <numeric/BodyPosition.fwd.hh>

// C++ headers
#include <utility/assert.hh>
#include <cmath>


namespace numeric {


/// @brief Unit quaternion 3-D orientation representation
///
/// @remarks
///  @li Quaternion defined as ( w, x, y, z )
///  @li Quaternions must be in the same coordinate frame to be used together
///  @li Successive rotation convention: q2 * q1 means q1 followed by q2
///  @li Performance-critical code can suppress normalization after every modification
///      and use explicit normalization
template< typename T >
class Quaternion
{


private: // Friends


	friend class BodyPosition< T >;


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


public: // Creation


	/// @brief Default constructor
	inline
	Quaternion() : // Sets to identity
		w_( T( 1 ) ),
		x_( T( 0 ) ),
		y_( T( 0 ) ),
		z_( T( 0 ) )
	{}


	/// @brief Coordinate constructor
	inline
	Quaternion(
		Value const & w_a,
		Value const & x_a,
		Value const & y_a,
		Value const & z_a,
		bool const precise = true
	) :
		w_( w_a ),
		x_( x_a ),
		y_( y_a ),
		z_( z_a )
	{
		assert( is_normalized() );
		if ( precise ) normalize();
	}


	/// @brief Identity named constructor
	inline
	static
	Quaternion
	identity()
	{
		return Quaternion( T( 1 ), T( 0 ), T( 0 ), T( 0 ) );
	}


	/// @brief Copy constructor
	inline
	Quaternion( Quaternion const & q ) :
		w_( q.w_ ),
		x_( q.x_ ),
		y_( q.y_ ),
		z_( q.z_ )
	{}


	/// @brief Destructor
	inline
	~Quaternion()
	{}


public: // copy assignment

	inline
	Quaternion &
	operator =( Quaternion const & q )
	{
		if ( this != &q ) {
			w_ = q.w_;
			x_ = q.x_;
			y_ = q.y_;
			z_ = q.z_;
		}
		return *this;
	}


public: // Properties


	/// @brief w
	inline
	Value const &
	w() const
	{
		return w_;
	}


	/// @brief x
	inline
	Value const &
	x() const
	{
		return x_;
	}


	/// @brief y
	inline
	Value const &
	y() const
	{
		return y_;
	}


	/// @brief z
	inline
	Value const &
	z() const
	{
		return z_;
	}


	/// @brief w squared
	inline
	Value
	w_squared() const
	{
		return w_ * w_;
	}


	/// @brief x squared
	inline
	Value
	x_squared() const
	{
		return x_ * x_;
	}


	/// @brief y squared
	inline
	Value
	y_squared() const
	{
		return y_ * y_;
	}


	/// @brief z squared
	inline
	Value
	z_squared() const
	{
		return z_ * z_;
	}


	/// @brief Norm: Should be one
	inline
	Value
	norm() const
	{
		return std::sqrt( ( w_ * w_ ) + ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Norm squared: Should be one
	inline
	Value
	norm_squared() const
	{
		return ( ( w_ * w_ ) + ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Norm error
	inline
	Value
	norm_error() const
	{
		return std::abs( T( 1 ) - norm() );
	}


	/// @brief Norm squared error
	inline
	Value
	norm_squared_error() const
	{
		return std::abs( T( 1 ) - norm_squared() );
	}


	/// @brief Magnitude: Should be one
	inline
	Value
	magnitude() const
	{
		return std::sqrt( ( w_ * w_ ) + ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Magnitude squared: Should be one
	inline
	Value
	magnitude_squared() const
	{
		return ( ( w_ * w_ ) + ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Magnitude error
	inline
	Value
	magnitude_error() const
	{
		return std::abs( T( 1 ) - magnitude() );
	}


	/// @brief Magnitude squared error
	inline
	Value
	magnitude_squared_error() const
	{
		return std::abs( T( 1 ) - magnitude_squared() );
	}


	/// @brief Magnitude squared error within tolerance?
	inline
	bool
	is_normalized( Value const & tol = Traits::quaternion_tolerance() ) const
	{
		return ( norm_squared_error() <= tol );
	}


	/// @brief Magnitude squared error exceeds tolerance?
	inline
	bool
	not_normalized( Value const & tol = Traits::quaternion_tolerance() ) const
	{
		return ( norm_squared_error() > tol );
	}


	/// @brief Principal angle of rotation (on [0,2*pi])
	inline
	Value
	angle() const
	{
		return T( 2 ) * std::acos( w_ );
	}


	/// @brief Axis of Rotation unit vector (direction for angle on [0,2*pi])
	inline
	Axis
	axis() const
	{
		return Axis( x_, y_, z_ ).normalize_or_zero(); // Returns zero vector if angle is zero
	}


	/// @brief Axis of rotation unit vector: Passed vector (direction for angle on [0,2*pi])
	inline
	Axis &
	axis( Axis & u ) const
	{
		return u.assign( x_, y_, z_ ).normalize_or_zero(); // Returns zero vector if angle is zero
	}


public: // Methods


	/// @brief Dot product
	inline
	Value
	dot( Quaternion const & q ) const
	{
		return ( w_ * q.w_ ) + ( x_ * q.x_ ) + ( y_ * q.y_ ) + ( z_ * q.z_ );
	}


	/// @brief Dot product
	inline
	Value
	dot_product( Quaternion const & q ) const
	{
		return ( w_ * q.w_ ) + ( x_ * q.x_ ) + ( y_ * q.y_ ) + ( z_ * q.z_ );
	}


public: // Methods: modifiers


	/// @brief Normalize
	inline
	Quaternion &
	normalize()
	{
		Value const norm_sq( norm_squared() );
		if ( norm_sq != T( 1 ) ) {
			assert( norm_sq > T( 0 ) );
			Value const norm_inv( T( 1 ) / std::sqrt( norm_sq ) );
			w_ *= norm_inv;
			x_ *= norm_inv;
			y_ *= norm_inv;
			z_ *= norm_inv;
		}
		return *this;
	}


	/// @brief Normalize if magnitude squared error exceeds tolerance
	inline
	Quaternion &
	normalize_if_needed( Value const & tol = Traits::quaternion_tolerance() )
	{
		Value const norm_sq( norm_squared() );
		if ( std::abs( T( 1 ) - norm_sq ) > tol ) {
			assert( norm_sq > T( 0 ) );
			Value const norm_inv( T( 1 ) / std::sqrt( norm_sq ) );
			w_ *= norm_inv;
			x_ *= norm_inv;
			y_ *= norm_inv;
			z_ *= norm_inv;
		}
		return *this;
	}


	/// @brief Identity
	inline
	Quaternion &
	to_identity()
	{
		w_ = T( 1 );
		x_ = T( 0 );
		y_ = T( 0 );
		z_ = T( 0 );
		return *this;
	}


	/// @brief Conjugate
	inline
	Quaternion &
	conjugate()
	{
		x_ = -x_;
		y_ = -y_;
		z_ = -z_;
		return *this;
	}


	/// @brief Invert
	inline
	Quaternion &
	invert()
	{
		x_ = -x_;
		y_ = -y_;
		z_ = -z_;
		return *this;
	}


	/// @brief Apply a successive Quaternion
	inline
	Quaternion &
	apply( Quaternion const & q, bool const precise = true )
	{
		return left_multiply_by( q, precise );
	}


	/// @brief Left multiply by a Quaternion
	inline
	Quaternion &
	left_multiply_by( Quaternion const & q, bool const precise = true )
	{
		Value const w_o( w_ );
		Value const x_o( x_ );
		Value const y_o( y_ );
		w_ = ( q.w_ * w_o ) - ( q.x_ * x_o ) - ( q.y_ * y_o ) - ( q.z_ * z_  );
		x_ = ( q.w_ * x_o ) + ( q.x_ * w_o ) + ( q.y_ * z_  ) - ( q.z_ * y_o );
		y_ = ( q.w_ * y_o ) - ( q.x_ * z_  ) + ( q.y_ * w_o ) + ( q.z_ * x_o );
		z_ = ( q.w_ * z_  ) + ( q.x_ * y_o ) - ( q.y_ * x_o ) + ( q.z_ * w_o );
		if ( precise ) normalize();
		return *this;
	}


	/// @brief Right multiply by a Quaternion
	inline
	Quaternion &
	right_multiply_by( Quaternion const & q, bool const precise = true )
	{
		Value const w_o( w_ );
		Value const x_o( x_ );
		Value const y_o( y_ );
		w_= ( w_o * q.w_ ) - ( x_o * q.x_ ) - ( y_o * q.y_ ) - ( z_ * q.z_ );
		x_= ( w_o * q.x_ ) + ( x_o * q.w_ ) + ( y_o * q.z_ ) - ( z_ * q.y_ );
		y_= ( w_o * q.y_ ) - ( x_o * q.z_ ) + ( y_o * q.w_ ) + ( z_ * q.x_ );
		z_= ( w_o * q.z_ ) + ( x_o * q.y_ ) - ( y_o * q.x_ ) + ( z_ * q.w_ );
		if ( precise ) normalize();
		return *this;
	}


	/// @brief Left multiply by the inverse of a Quaternion
	inline
	Quaternion &
	left_multiply_by_inverse_of( Quaternion const & q, bool const precise = true )
	{
		Value const w_o( w_ );
		Value const x_o( x_ );
		Value const y_o( y_ );
		w_ = ( q.w_ * w_o ) - ( q.x_ * x_o ) - ( q.y_ * y_o ) - ( q.z_ * z_  );
		x_ = ( q.w_ * x_o ) + ( q.x_ * w_o ) + ( q.y_ * z_  ) - ( q.z_ * y_o );
		y_ = ( q.w_ * y_o ) - ( q.x_ * z_  ) + ( q.y_ * w_o ) + ( q.z_ * x_o );
		z_ = ( q.w_ * z_  ) + ( q.x_ * y_o ) - ( q.y_ * x_o ) + ( q.z_ * w_o );
		if ( precise ) normalize();
		return *this;
	}


	/// @brief Right multiply by the inverse of a Quaternion
	inline
	Quaternion &
	right_multiply_by_inverse_of( Quaternion const & q, bool const precise = true )
	{
		Value const w_o( w_ );
		Value const x_o( x_ );
		Value const y_o( y_ );
		w_ = ( w_o * q.w_ ) - ( x_o * q.x_ ) - ( y_o * q.y_ ) - ( z_ * q.z_ );
		x_ = ( w_o * q.x_ ) + ( x_o * q.w_ ) + ( y_o * q.z_ ) - ( z_ * q.y_ );
		y_ = ( w_o * q.y_ ) - ( x_o * q.z_ ) + ( y_o * q.w_ ) + ( z_ * q.x_ );
		z_ = ( w_o * q.z_ ) + ( x_o * q.y_ ) - ( y_o * q.x_ ) + ( z_ * q.w_ );
		if ( precise ) normalize();
		return *this;
	}


	/// @brief Swap
	inline
	void
	swap( Quaternion & q )
	{
		Value t;
		t = w_; w_ = q.w_; q.w_ = t;
		t = x_; x_ = q.x_; q.x_ = t;
		t = y_; y_ = q.y_; q.y_ = t;
		t = z_; z_ = q.z_; q.z_ = t;
	}


public: // Methods: generators


	/// @brief Conjugated
	inline
	Quaternion
	conjugated() const
	{
		return Quaternion( w_, -x_, -y_, -z_ );
	}


	/// @brief Inverse
	inline
	Quaternion
	inverse() const
	{
		return Quaternion( w_, -x_, -y_, -z_ );
	}


	/// @brief Quaternion * Quaternion
	friend
	inline
	Quaternion
	operator *( Quaternion const & q2, Quaternion const & q1 )
	{
		return Quaternion(
			( q2.w_ * q1.w_ ) - ( q2.x_ * q1.x_ ) - ( q2.y_ * q1.y_ ) - ( q2.z_ * q1.z_ ),
			( q2.w_ * q1.x_ ) + ( q2.x_ * q1.w_ ) + ( q2.y_ * q1.z_ ) - ( q2.z_ * q1.y_ ),
			( q2.w_ * q1.y_ ) - ( q2.x_ * q1.z_ ) + ( q2.y_ * q1.w_ ) + ( q2.z_ * q1.x_ ),
			( q2.w_ * q1.z_ ) + ( q2.x_ * q1.y_ ) - ( q2.y_ * q1.x_ ) + ( q2.z_ * q1.w_ )
			).normalize();
	}


	/// @brief Product: Quaternion * Quaternion
	friend
	inline
	Quaternion
	product( Quaternion const & q2, Quaternion const & q1, bool const precise = true )
	{
		if ( precise ) { // Return normalized product
			return Quaternion(
				( q2.w_ * q1.w_ ) - ( q2.x_ * q1.x_ ) - ( q2.y_ * q1.y_ ) - ( q2.z_ * q1.z_ ),
				( q2.w_ * q1.x_ ) + ( q2.x_ * q1.w_ ) + ( q2.y_ * q1.z_ ) - ( q2.z_ * q1.y_ ),
				( q2.w_ * q1.y_ ) - ( q2.x_ * q1.z_ ) + ( q2.y_ * q1.w_ ) + ( q2.z_ * q1.x_ ),
				( q2.w_ * q1.z_ ) + ( q2.x_ * q1.y_ ) - ( q2.y_ * q1.x_ ) + ( q2.z_ * q1.w_ )
				).normalize();
		} else { // Return product
			return Quaternion(
				( q2.w_ * q1.w_ ) - ( q2.x_ * q1.x_ ) - ( q2.y_ * q1.y_ ) - ( q2.z_ * q1.z_ ),
				( q2.w_ * q1.x_ ) + ( q2.x_ * q1.w_ ) + ( q2.y_ * q1.z_ ) - ( q2.z_ * q1.y_ ),
				( q2.w_ * q1.y_ ) - ( q2.x_ * q1.z_ ) + ( q2.y_ * q1.w_ ) + ( q2.z_ * q1.x_ ),
				( q2.w_ * q1.z_ ) + ( q2.x_ * q1.y_ ) - ( q2.y_ * q1.x_ ) + ( q2.z_ * q1.w_ )
			);
		}
	}


	/// @brief Identity Quaternion for expressions
	/// @note  Default and identity() named constructors can be faster for construction
	/// @note  Can be safely used in construction of global objects
	inline
	static
	Quaternion const &
	I()
	{
		static Quaternion const I_( T( 1 ), T( 0 ), T( 0 ), T( 0 ) );
		return I_;
	}


public: // Comparison


	/// @brief Quaternion == Quaternion
	friend
	inline
	bool
	operator ==( Quaternion const & q1, Quaternion const & q2 )
	{
		return ( ( q1.w_ == q2.w_ ) && ( q1.x_ == q2.x_ ) && ( q1.y_ == q2.y_ ) && ( q1.z_ == q2.z_ ) );
	}


	/// @brief Quaternion != Quaternion
	friend
	inline
	bool
	operator !=( Quaternion const & q1, Quaternion const & q2 )
	{
		return ( ( q1.w_ != q2.w_ ) || ( q1.x_ != q2.x_ ) || ( q1.y_ != q2.y_ ) || ( q1.z_ != q2.z_ ) );
	}


	/// @brief Dot product
	friend
	inline
	Value
	dot( Quaternion const & q1, Quaternion const & q2 )
	{
		return ( q1.w_ * q2.w_ ) + ( q1.x_ * q2.x_ ) + ( q1.y_ * q2.y_ ) + ( q1.z_ * q2.z_ );
	}


	/// @brief Dot product
	friend
	inline
	Value
	dot_product( Quaternion const & q1, Quaternion const & q2 )
	{
		return ( q1.w_ * q2.w_ ) + ( q1.x_ * q2.x_ ) + ( q1.y_ * q2.y_ ) + ( q1.z_ * q2.z_ );
	}


private: // Fields


	/// @brief w coordinate
	Value w_;

	/// @brief x coordinate
	Value x_;

	/// @brief y coordinate
	Value y_;

	/// @brief z coordinate
	Value z_;


}; // Quaternion


/// @brief Quaternion * Quaternion
template< typename T >
Quaternion< T >
operator *( Quaternion< T > const & q2, Quaternion< T > const & q1 );


/// @brief Product: Quaternion * Quaternion
template< typename T >
Quaternion< T >
product( Quaternion< T > const & q2, Quaternion< T > const & q1, bool const precise );


/// @brief Quaternion == Quaternion
template< typename T >
bool
operator ==( Quaternion< T > const & q1, Quaternion< T > const & q2 );


/// @brief Quaternion != Quaternion
template< typename T >
bool
operator !=( Quaternion< T > const & q1, Quaternion< T > const & q2 );


/// @brief Dot product
template< typename T >
T
dot( Quaternion< T > const & q1, Quaternion< T > const & q2 );


/// @brief Dot product
template< typename T >
T
dot_product( Quaternion< T > const & q1, Quaternion< T > const & q2 );


} // namespace numeric


#endif // INCLUDED_numeric_Quaternion_HH
