// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzVector.hh
/// @brief  Fast (x,y,z)-coordinate numeric vector
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Inline, loop-free functions for speed
///  @li Non-virtual destructor for speed: Not set up for use as a base class
///  @li Pointer constructor and assignment not available for xyzVectors of pointers
///  @li Numeric vector semantics: spatial partial ordering


#ifndef INCLUDED_numeric_xyzVector_hh
#define INCLUDED_numeric_xyzVector_hh


// Unit headers
#include <numeric/xyzVector.fwd.hh>

// Package headers
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/trig.functions.hh>
#include <platform/types.hh>

// C++ headers
#include <stdexcept>
#include <utility/assert.hh>
#include <cmath>
#ifdef GL_GRAPHICS
	#include <numeric/xyz.functions.hh>
#endif

#include <ObjexxFCL/FArrayTraits.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Boost Headers
#include <boost/functional/hash.hpp>

// C++ headers
#include <cstdlib>

namespace numeric {


/// @brief xyzVector: Fast (x,y,z)-coordinate numeric vector
template< typename T >
class xyzVector
{


private: // Friends


	template< typename > friend class xyzVector;
	template< typename > friend class xyzMatrix;

	// Friend functions (for speed of non-inlining debug builds)
	friend xyzVector< T > operator *<>( xyzMatrix< T > const & m, xyzVector< T > const & v );
	friend xyzVector< T > product<>( xyzMatrix< T > const & m, xyzVector< T > const & v );
	friend xyzVector< T > & inplace_product<>( xyzMatrix< T > const & m, xyzVector< T > & v );
	friend xyzVector< T > transpose_product<>( xyzMatrix< T > const & m, xyzVector< T > const & v );
	friend xyzVector< T > & inplace_transpose_product<>( xyzMatrix< T > const & m, xyzVector< T > & v );
	friend xyzMatrix< T > outer_product<>( xyzVector< T > const & a, xyzVector< T > const & b );
	friend xyzMatrix< T > projection_matrix<>( xyzVector< T > const & v );
	friend xyzMatrix< T > rotation_matrix<>( xyzVector< T > const & axis, T const & theta );


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

	// Types to prevent compile failure when std::distance is in scope
	typedef  void  iterator_category;
	typedef  void  difference_type;


public: // Creation


	/// @brief Default constructor
	/// @note  Values are uninitialized for efficiency
	inline
	xyzVector()
	{}


	/// @brief Copy constructor
	inline
	xyzVector( xyzVector const & v ) :
		x_( v.x_ ),
		y_( v.y_ ),
		z_( v.z_ )
	{}


	/// @brief Copy constructor
	template< typename U >
	inline
	xyzVector( xyzVector< U > const & v ) :
		x_( v.x_ ),
		y_( v.y_ ),
		z_( v.z_ )
	{}


	/// @brief Uniform value constructor
	inline
	explicit
	xyzVector( Value const & t ) :
		x_( t ),
		y_( t ),
		z_( t )
	{}


	/// @brief Triple value constructor
	inline
	xyzVector(
		Value const & x_a,
		Value const & y_a,
		Value const & z_a
	 ) :
		x_( x_a ),
		y_( y_a ),
		z_( z_a )
	{}


	/// @brief Pointer to contiguous values constructor
	/// @note  U must be assignable to a Value
	/// @warning No way to check that argument points to three values
	/// @warning Argument missing an & operator will quietly call the uniform value constructor
	template< typename U >
	inline
	explicit
	xyzVector( U const * p ) :
		x_( *p ),
		y_( *++p ),
		z_( *++p )
	{}


	/// @brief Destructor
	inline
	~xyzVector()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	xyzVector &
	operator =( xyzVector const & v )
	{
		if ( this != &v ) {
			x_ = v.x_;
			y_ = v.y_;
			z_ = v.z_;
		}
		return *this;
	}


	/// @brief Copy assignment
	template< typename U >
	inline
	xyzVector &
	operator =( xyzVector< U > const & v )
	{
		x_ = v.x_;
		y_ = v.y_;
		z_ = v.z_;
		return *this;
	}


	/// @brief Assignment from pointer to contiguous values
	/// @warning No way to check that argument points to three values
	template< typename U >
	inline
	xyzVector &
	operator =( U const * p )
	{
		x_ = *p;
		y_ = *++p;
		z_ = *++p;
		return *this;
	}


	/// @brief += xyzVector
	template< typename U >
	inline
	xyzVector &
	operator +=( xyzVector< U > const & v )
	{
		x_ += v.x_;
		y_ += v.y_;
		z_ += v.z_;
		return *this;
	}


	/// @brief -= xyzVector
	template< typename U >
	inline
	xyzVector &
	operator -=( xyzVector< U > const & v )
	{
		x_ -= v.x_;
		y_ -= v.y_;
		z_ -= v.z_;
		return *this;
	}


	/// @brief Assign Value * xyzVector
	/// @note  Avoids temporary of = Value * xyzVector
	template< typename U >
	inline
	xyzVector &
	scaled_assign( Value const & t, xyzVector< U > const & v )
	{
		x_ = t * v.x_;
		y_ = t * v.y_;
		z_ = t * v.z_;
		return *this;
	}


	/// @brief Add Value * xyzVector
	/// @note  Avoids temporary of += Value * xyzVector
	template< typename U >
	inline
	xyzVector &
	scaled_add( Value const & t, xyzVector< U > const & v )
	{
		x_ += t * v.x_;
		y_ += t * v.y_;
		z_ += t * v.z_;
		return *this;
	}


	/// @brief Subtract Value * xyzVector
	/// @note  Avoids temporary of -= Value * xyzVector
	template< typename U >
	inline
	xyzVector &
	scaled_sub( Value const & t, xyzVector< U > const & v )
	{
		x_ -= t * v.x_;
		y_ -= t * v.y_;
		z_ -= t * v.z_;
		return *this;
	}


	/// @brief = Value
	inline
	xyzVector &
	operator =( Value const & t )
	{
		x_ = y_ = z_ = t;
		return *this;
	}


	/// @brief += Value
	inline
	xyzVector &
	operator +=( Value const & t )
	{
		x_ += t;
		y_ += t;
		z_ += t;
		return *this;
	}


	/// @brief -= Value
	inline
	xyzVector &
	operator -=( Value const & t )
	{
		x_ -= t;
		y_ -= t;
		z_ -= t;
		return *this;
	}


	/// @brief *= Value
	inline
	xyzVector &
	operator *=( Value const & t )
	{
		x_ *= t;
		y_ *= t;
		z_ *= t;
		return *this;
	}


	/// @brief /= Value
	inline
	xyzVector &
	operator /=( Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		x_ *= inv_t;
		y_ *= inv_t;
		z_ *= inv_t;
		return *this;
	}


	/// @brief Triple value assignment
	inline
	xyzVector &
	assign(
		Value const & x_a,
		Value const & y_a,
		Value const & z_a
	)
	{
		x_ = x_a;
		y_ = y_a;
		z_ = z_a;
		return *this;
	}


public: // Methods
	
    
    /// @brief to_string, useful for utility exits
	inline
	std::string
	to_string() const
	{
		return "(" + utility::to_string(x_) + ", " + utility::to_string(y_) + ", " + utility::to_string(z_) + ")";
	}

	
	/// @brief Show
	inline
	void
	show(std::ostream & output=std::cout) const
	{
		output << "(" << x_ << ", " << y_ << ", " << z_ << ")";
	}


	/// @brief Clear
	inline
	xyzVector &
	clear()
	{
		x_ = y_ = z_ = Value( 0 );
		return *this;
	}


	/// @brief Zero
	inline
	xyzVector &
	zero()
	{
		x_ = y_ = z_ = Value( 0 );
		return *this;
	}


	/// @brief Negate
	inline
	xyzVector &
	negate()
	{
		x_ = -x_;
		y_ = -y_;
		z_ = -z_;
		return *this;
	}


	/// @brief -xyzVector (negated copy)
	inline
	xyzVector
	operator -() const
	{
		return xyzVector( -x_, -y_, -z_ );
	}


	/// @brief Negated copy
	inline
	xyzVector
	negated() const
	{
		return xyzVector( -x_, -y_, -z_ );
	}


	/// @brief Negated: Return via argument (slightly faster)
	inline
	void
	negated( xyzVector & a ) const
	{
		a.x_ = -x_;
		a.y_ = -y_;
		a.z_ = -z_;
	}


	/// @brief xyzVector + xyzVector
	friend
	inline
	xyzVector
	operator +( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector( a.x_ + b.x_, a.y_ + b.y_, a.z_ + b.z_ );
	}


	/// @brief xyzVector + Value
	friend
	inline
	xyzVector
	operator +( xyzVector const & v, Value const & t )
	{
		return xyzVector( v.x_ + t, v.y_ + t, v.z_ + t );
	}


	/// @brief Value + xyzVector
	friend
	inline
	xyzVector
	operator +( Value const & t, xyzVector const & v )
	{
		return xyzVector( t + v.x_, t + v.y_, t + v.z_ );
	}


	/// @brief xyzVector - xyzVector
	friend
	inline
	xyzVector
	operator -( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector( a.x_ - b.x_, a.y_ - b.y_, a.z_ - b.z_ );
	}


	/// @brief xyzVector - Value
	friend
	inline
	xyzVector
	operator -( xyzVector const & v, Value const & t )
	{
		return xyzVector( v.x_ - t, v.y_ - t, v.z_ - t );
	}


	/// @brief Value - xyzVector
	friend
	inline
	xyzVector
	operator -( Value const & t, xyzVector const & v )
	{
		return xyzVector( t - v.x_, t - v.y_, t - v.z_ );
	}


	/// @brief xyzVector * Value
	friend
	inline
	xyzVector
	operator *( xyzVector const & v, Value const & t )
	{
		return xyzVector( v.x_ * t, v.y_ * t, v.z_ * t );
	}


	/// @brief Value * xyzVector
	friend
	inline
	xyzVector
	operator *( Value const & t, xyzVector const & v )
	{
		return xyzVector( t * v.x_, t * v.y_, t * v.z_ );
	}


	/// @brief xyzVector / Value
	friend
	inline
	xyzVector
	operator /( xyzVector const & v, Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value ( 1 ) / t );
		return xyzVector( v.x_ * inv_t, v.y_ * inv_t, v.z_ * inv_t );
	}


	/// @brief Add: xyzVector + xyzVector
	friend
	inline
	void
	add( xyzVector const & a, xyzVector const & b, xyzVector & r )
	{
		r.x_ = a.x_ + b.x_;
		r.y_ = a.y_ + b.y_;
		r.z_ = a.z_ + b.z_;
	}


	/// @brief Add: xyzVector + Value
	friend
	inline
	void
	add( xyzVector const & v, Value const & t, xyzVector & r )
	{
		r.x_ = v.x_ + t;
		r.y_ = v.y_ + t;
		r.z_ = v.z_ + t;
	}


	/// @brief Add: Value + xyzVector
	friend
	inline
	void
	add( Value const & t, xyzVector const & v, xyzVector & r )
	{
		r.x_ = t + v.x_;
		r.y_ = t + v.y_;
		r.z_ = t + v.z_;
	}


	/// @brief Subtract: xyzVector - xyzVector
	friend
	inline
	void
	subtract( xyzVector const & a, xyzVector const & b, xyzVector & r )
	{
		r.x_ = a.x_ - b.x_;
		r.y_ = a.y_ - b.y_;
		r.z_ = a.z_ - b.z_;
	}


	/// @brief Subtract: xyzVector - Value
	friend
	inline
	void
	subtract( xyzVector const & v, Value const & t, xyzVector & r )
	{
		r.x_ = v.x_ - t;
		r.y_ = v.y_ - t;
		r.z_ = v.z_ - t;
	}


	/// @brief Subtract: Value - xyzVector
	friend
	inline
	void
	subtract( Value const & t, xyzVector const & v, xyzVector & r )
	{
		r.x_ = t - v.x_;
		r.y_ = t - v.y_;
		r.z_ = t - v.z_;
	}


	/// @brief Multiply: xyzVector * Value
	friend
	inline
	void
	multiply( xyzVector const & v, Value const & t, xyzVector & r )
	{
		r.x_ = v.x_ * t;
		r.y_ = v.y_ * t;
		r.z_ = v.z_ * t;
	}


	/// @brief Multiply: Value * xyzVector
	friend
	inline
	void
	multiply( Value const & t, xyzVector const & v, xyzVector & r )
	{
		r.x_ = t * v.x_;
		r.y_ = t * v.y_;
		r.z_ = t * v.z_;
	}


	/// @brief Divide: xyzVector / Value
	friend
	inline
	void
	divide( xyzVector const & v, Value const & t, xyzVector & r )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		r.x_ = v.x_ * inv_t;
		r.y_ = v.y_ * inv_t;
		r.z_ = v.z_ * inv_t;
	}


	/// @brief Set minimum coordinates wrt another xyzVector
	inline
	xyzVector &
	min( xyzVector const & v )
	{
		x_ = ( x_ <= v.x_ ? x_ : v.x_ );
		y_ = ( y_ <= v.y_ ? y_ : v.y_ );
		z_ = ( z_ <= v.z_ ? z_ : v.z_ );
		return *this;
	}


	/// @brief Set maximum coordinates wrt another xyzVector
	inline
	xyzVector &
	max( xyzVector const & v )
	{
		x_ = ( x_ >= v.x_ ? x_ : v.x_ );
		y_ = ( y_ >= v.y_ ? y_ : v.y_ );
		z_ = ( z_ >= v.z_ ? z_ : v.z_ );
		return *this;
	}


	/// @brief xyzVector with min coordinates of two xyzVectors
	friend
	inline
	xyzVector
	min( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
		 ( a.x_ <= b.x_ ? a.x_ : b.x_ ),
		 ( a.y_ <= b.y_ ? a.y_ : b.y_ ),
		 ( a.z_ <= b.z_ ? a.z_ : b.z_ )
		);
	}


	/// @brief xyzVector with max coordinates of two xyzVectors
	friend
	inline
	xyzVector
	max( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
		 ( a.x_ >= b.x_ ? a.x_ : b.x_ ),
		 ( a.y_ >= b.y_ ? a.y_ : b.y_ ),
		 ( a.z_ >= b.z_ ? a.z_ : b.z_ )
		);
	}


	/// @brief Normalize
	inline
	xyzVector &
	normalize()
	{
		Value const length_ = length();
		//assert ( length_ != Value ( 0 ));
		if ( length_ == Value ( 0 )){
			throw utility::excn::EXCN_BadInput("Cannot normalize xyzVector of length() zero");
		}
		
		Value const inv_length( Value( 1 ) / length_ );
		x_ *= inv_length;
		y_ *= inv_length;
		z_ *= inv_length;
		return *this;
	}

	/// @brief get the lowest of the three values in the vector
	inline
	T
	minimum_value() const
	{
		return std::min(std::min(x_,y_),std::min(y_,z_));
	}

	/// @brief get the highest of the three values in the vector
	inline
	T
	maximum_value() const
	{
		return std::max(std::max(x_,y_),std::max(y_,z_));
	}


	/// @brief Normalized
	inline
	void
	normalized( xyzVector & a ) const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const inv_length( Value( 1 ) / length_ );
		a.x_ = x_ * inv_length;
		a.y_ = y_ * inv_length;
		a.z_ = z_ * inv_length;
	}


	/// @brief Normalize: zero xyzVector if length is zero
	inline
	xyzVector &
	normalize_or_zero()
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			x_ *= inv_length;
			y_ *= inv_length;
			z_ *= inv_length;
		} else { // Set zero vector
			x_ = Value( 0 );
			y_ = Value( 0 );
			z_ = Value( 0 );
		}
		return *this;
	}


	/// @brief Normalized: zero xyzVector if length is zero
	inline
	void
	normalized_or_zero( xyzVector & a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			a.x_ = x_ * inv_length;
			a.y_ = y_ * inv_length;
			a.z_ = z_ * inv_length;
		} else { // Return zero vector
			a.x_ = Value( 0 );
			a.y_ = Value( 0 );
			a.z_ = Value( 0 );
		}
	}


	/// @brief Normalize: arbitrary normalized xyzVector if length is zero
	inline
	xyzVector &
	normalize_any()
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			x_ *= inv_length;
			y_ *= inv_length;
			z_ *= inv_length;
		} else { // Set arbitrary normalized vector
			x_ = Value( 1 );
			y_ = Value( 0 );
			z_ = Value( 0 );
		}
		return *this;
	}


	/// @brief Normalized: arbitrary normalized xyzVector if length is zero
	inline
	void
	normalized_any( xyzVector & a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			a.x_ = x_ * inv_length;
			a.y_ = y_ * inv_length;
			a.z_ = z_ * inv_length;
		} else { // Return arbitrary normalized vector
			a.x_ = Value( 1 );
			a.y_ = Value( 0 );
			a.z_ = Value( 0 );
		}
	}


	/// @brief Normalize to a length
	inline
	xyzVector &
	normalize( Value const & length_a )
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const dilation = length_a / length_;
		x_ *= dilation;
		y_ *= dilation;
		z_ *= dilation;
		return *this;
	}


	/// @brief Normalized to a length
	inline
	void
	normalized( Value const & length_a, xyzVector & a ) const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const dilation = length_a / length_;
		a.x_ = x_ * dilation;
		a.y_ = y_ * dilation;
		a.z_ = z_ * dilation;
	}


	/// @brief Normalize to a length: zero xyzVector if length is zero
	inline
	xyzVector &
	normalize_or_zero( Value const & length_a )
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			x_ *= dilation;
			y_ *= dilation;
			z_ *= dilation;
		} else { // Set zero vector
			x_ = Value( 0 );
			y_ = Value( 0 );
			z_ = Value( 0 );
		}
		return *this;
	}


	/// @brief Normalized to a length: zero xyzVector if length is zero
	inline
	void
	normalized_or_zero( Value const & length_a, xyzVector & a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			a.x_ = x_ * dilation;
			a.y_ = y_ * dilation;
			a.z_ = z_ * dilation;
		} else { // Return zero vector
			a.x_ = Value( 0 );
			a.y_ = Value( 0 );
			a.z_ = Value( 0 );
		}
	}


	/// @brief Normalize to a length: arbitrary normalized xyzVector if length is zero
	inline
	xyzVector &
	normalize_any( Value const & length_a )
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			x_ *= dilation;
			y_ *= dilation;
			z_ *= dilation;
		} else { // Set arbitrary normalized vector
			x_ = length_a;
			y_ = Value( 0 );
			z_ = Value( 0 );
		}
		return *this;
	}


	/// @brief Normalized to a length: arbitrary normalized xyzVector if length is zero
	inline
	void
	normalized_any( Value const & length_a, xyzVector & a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			a.x_ = x_ * dilation;
			a.y_ = y_ * dilation;
			a.z_ = z_ * dilation;
		} else { // Return arbitrary normalized vector
			a.x_ = length_a;
			a.y_ = Value( 0 );
			a.z_ = Value( 0 );
		}
	}


	/// @brief Normalized copy
	inline
	xyzVector
	normalized() const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const inv_length( Value( 1 ) / length_ );
		return xyzVector( x_ * inv_length, y_ * inv_length, z_ * inv_length );
	}


	/// @brief Normalized copy: Zero xyzVector if length is zero
	inline
	xyzVector
	normalized_or_zero() const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			return xyzVector( x_ * inv_length, y_ * inv_length, z_ * inv_length );
		} else { // Return zero vector
			return xyzVector( Value( 0 ), Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Normalized copy: Arbitrary normalized xyzVector if length is zero
	inline
	xyzVector
	normalized_any() const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			return xyzVector( x_ * inv_length, y_ * inv_length, z_ * inv_length );
		} else { // Return arbitrary normalized vector
			return xyzVector( Value( 1 ), Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Normalized to a length copy
	inline
	xyzVector
	normalized( Value const & length_a ) const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const dilation = length_a / length_;
		return xyzVector( x_ * dilation, y_ * dilation, z_ * dilation );
	}


	/// @brief Normalized to a length copy: Zero xyzVector if length is zero
	inline
	xyzVector
	normalized_or_zero( Value const & length_a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			return xyzVector( x_ * dilation, y_ * dilation, z_ * dilation );
		} else { // Return zero vector
			return xyzVector( Value( 0 ), Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Normalized to a length copy: Arbitrary normalized xyzVector if length is zero
	inline
	xyzVector
	normalized_any( Value const & length_a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			return xyzVector( x_ * dilation, y_ * dilation, z_ * dilation );
		} else { // Return arbitrary normalized vector
			return xyzVector( length_a, Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Project normal
	/// @note  This vector projected normally onto input vector
	/// @note  Not meaningful when v == 0
	inline
	xyzVector &
	project_normal( xyzVector const & v )
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		x_ -= c * v.x_;
		y_ -= c * v.y_;
		z_ -= c * v.z_;
		return *this;
	}


	/// @brief Projected normal copy
	/// @note  Copy of this vector projected normally onto input vector
	/// @note  Not meaningful when v == 0
	inline
	xyzVector
	projected_normal( xyzVector const & v ) const
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		return xyzVector( x_ - ( c * v.x_ ), y_ - ( c * v.y_ ), z_ - ( c * v.z_ ) );
	}


	/// @brief Projected normal
	/// @note  Copy of this vector projected normally onto first input vector
	/// @note  Not meaningful when v == 0
	inline
	void
	projected_normal( xyzVector const & v, xyzVector & a ) const
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		a.x_ = x_ - ( c * v.x_ );
		a.y_ = y_ - ( c * v.y_ );
		a.z_ = z_ - ( c * v.z_ );
	}


	/// @brief Project parallel
	/// @note  This vector projected in direction of input vector
	/// @note  Not meaningful when v == 0
	inline
	xyzVector &
	project_parallel( xyzVector const & v )
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		x_ = c * v.x_;
		y_ = c * v.y_;
		z_ = c * v.z_;
		return *this;
	}


	/// @brief Projected parallel copy
	/// @note  Copy of this vector projected in direction of input vector
	/// @note  Not meaningful when v == 0
	inline
	xyzVector
	projected_parallel( xyzVector const & v ) const
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		return xyzVector( c * v.x_, c * v.y_, c * v.z_ );
	}


	/// @brief Projected parallel
	/// @note  Copy of this vector projected in direction of first input vector
	/// @note  Not meaningful when v == 0
	inline
	void
	projected_parallel( xyzVector const & v, xyzVector & a )
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		a.x_ = c * v.x_;
		a.y_ = c * v.y_;
		a.z_ = c * v.z_;
	}


	/// @brief Distance
	inline
	Value
	distance( xyzVector const & v ) const
	{
		return std::sqrt( square( x_ - v.x_ ) + square( y_ - v.y_ ) + square( z_ - v.z_ ) );
	}


	/// @brief Distance
//Commented out 8/9/11 SML; xyzVector friend functions cause problems for some versions of GCC.  Read the mailing list logs for 2/17/09 for further information.
	// friend
	// inline
	// Value
	// distance( xyzVector const & a, xyzVector const & b )
	// {
	// 	return std::sqrt( square( a.x_ - b.x_ ) + square( a.y_ - b.y_ ) + square( a.z_ - b.z_ ) );
	// }


	/// @brief Distance squared
	inline
	Value
	distance_squared( xyzVector const & v ) const
	{
		return square( x_ - v.x_ ) + square( y_ - v.y_ ) + square( z_ - v.z_ );
	}


	/// @brief Distance squared
//Commented out 8/9/11 SML; xyzVector friend functions cause problems for some versions of GCC.  Read the mailing list logs for 2/17/09 for further information.
	// friend
	// inline
	// Value
	// distance_squared( xyzVector const & a, xyzVector const & b )
	// {
	// 	return square( a.x_ - b.x_ ) + square( a.y_ - b.y_ ) + square( a.z_ - b.z_ );
	// }


	/// @brief Dot product
	inline
	Value
	dot( xyzVector const & v ) const
	{
		return ( x_ * v.x_ ) + ( y_ * v.y_ ) + ( z_ * v.z_ );
	}


	/// @brief Dot product
	inline
	Value
	dot_product( xyzVector const & v ) const
	{
		return ( x_ * v.x_ ) + ( y_ * v.y_ ) + ( z_ * v.z_ );
	}


	/// @brief Inner product ( == dot product )
	inline
	Value
	inner_product( xyzVector const & v ) const
	{
		return ( x_ * v.x_ ) + ( y_ * v.y_ ) + ( z_ * v.z_ );
	}


	/// @brief Dot product
	friend
	inline
	Value
	dot( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ * b.x_ ) + ( a.y_ * b.y_ ) + ( a.z_ * b.z_ );
	}


	/// @brief Dot product
	friend
	inline
	Value
	dot_product( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ * b.x_ ) + ( a.y_ * b.y_ ) + ( a.z_ * b.z_ );
	}


	/// @brief Inner product ( == dot product )
	friend
	inline
	Value
	inner_product( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ * b.x_ ) + ( a.y_ * b.y_ ) + ( a.z_ * b.z_ );
	}


	/// @brief Cross product
	inline
	xyzVector
	cross( xyzVector const & v ) const
	{
		return xyzVector(
			( y_ * v.z_ ) - ( z_ * v.y_ ),
			( z_ * v.x_ ) - ( x_ * v.z_ ),
			( x_ * v.y_ ) - ( y_ * v.x_ )
		);
	}


	/// @brief Cross product
	inline
	xyzVector
	cross_product( xyzVector const & v ) const
	{
		return xyzVector(
			( y_ * v.z_ ) - ( z_ * v.y_ ),
			( z_ * v.x_ ) - ( x_ * v.z_ ),
			( x_ * v.y_ ) - ( y_ * v.x_ )
		);
	}
	
	/// @brief Cross product
	friend
	inline
	xyzVector
	cross( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
			( a.y_ * b.z_ ) - ( a.z_ * b.y_ ),
			( a.z_ * b.x_ ) - ( a.x_ * b.z_ ),
			( a.x_ * b.y_ ) - ( a.y_ * b.x_ )
		);
	}


	/// @brief Cross product
	friend
	inline
	xyzVector
	cross_product( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
			( a.y_ * b.z_ ) - ( a.z_ * b.y_ ),
			( a.z_ * b.x_ ) - ( a.x_ * b.z_ ),
			( a.x_ * b.y_ ) - ( a.y_ * b.x_ )
		);
	}


	/// @brief Cross product: Return via argument (slightly faster)
	friend
	inline
	void
	cross( xyzVector const & a, xyzVector const & b, xyzVector & c )
	{
		c.x_ = ( a.y_ * b.z_ ) - ( a.z_ * b.y_ );
		c.y_ = ( a.z_ * b.x_ ) - ( a.x_ * b.z_ );
		c.z_ = ( a.x_ * b.y_ ) - ( a.y_ * b.x_ );
	}


	/// @brief Cross product: Return via argument (slightly faster)
	friend
	inline
	void
	cross_product( xyzVector const & a, xyzVector const & b, xyzVector & c )
	{
		c.x_ = ( a.y_ * b.z_ ) - ( a.z_ * b.y_ );
		c.y_ = ( a.z_ * b.x_ ) - ( a.x_ * b.z_ );
		c.z_ = ( a.x_ * b.y_ ) - ( a.y_ * b.x_ );
	}


	/// @brief Midpoint of 2 xyzVectors
	friend
	inline
	xyzVector
	midpoint( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
			Value( 0.5 * ( a.x_ + b.x_ ) ),
			Value( 0.5 * ( a.y_ + b.y_ ) ),
			Value( 0.5 * ( a.z_ + b.z_ ) )
		);
	}


	/// @brief Midpoint of 2 xyzVectors: Return via argument (slightly faster)
	friend
	inline
	void
	midpoint( xyzVector const & a, xyzVector const & b, xyzVector & m )
	{
		m.x_ = Value( 0.5 * ( a.x_ + b.x_ ) );
		m.y_ = Value( 0.5 * ( a.y_ + b.y_ ) );
		m.z_ = Value( 0.5 * ( a.z_ + b.z_ ) );
	}


	/// @brief Center of 2 xyzVectors
	friend
	inline
	xyzVector
	center( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
			Value( 0.5 * ( a.x_ + b.x_ ) ),
			Value( 0.5 * ( a.y_ + b.y_ ) ),
			Value( 0.5 * ( a.z_ + b.z_ ) )
		);
	}


	/// @brief Center of 2 xyzVectors: Return via argument (slightly faster)
	friend
	inline
	void
	center( xyzVector const & a, xyzVector const & b, xyzVector & m )
	{
		m.x_ = Value( 0.5 * ( a.x_ + b.x_ ) );
		m.y_ = Value( 0.5 * ( a.y_ + b.y_ ) );
		m.z_ = Value( 0.5 * ( a.z_ + b.z_ ) );
	}


	/// @brief Center of 3 xyzVectors
	friend
	inline
	xyzVector
	center( xyzVector const & a, xyzVector const & b, xyzVector const & c )
	{
		long double const third( 1.0 / 3.0 );
		return xyzVector(
			Value( third * ( a.x_ + b.x_ + c.x_ ) ),
			Value( third * ( a.y_ + b.y_ + c.y_ ) ),
			Value( third * ( a.z_ + b.z_ + c.z_ ) )
		);
	}


	/// @brief Center of 3 xyzVectors: Return via argument (slightly faster)
	friend
	inline
	void
	center( xyzVector const & a, xyzVector const & b, xyzVector const & c, xyzVector & m )
	{
		long double const third( 1.0 / 3.0 );
		m.x_ = Value( third * ( a.x_ + b.x_ + c.x_ ) );
		m.y_ = Value( third * ( a.y_ + b.y_ + c.y_ ) );
		m.z_ = Value( third * ( a.z_ + b.z_ + c.z_ ) );
	}


	/// @brief Center of 4 xyzVectors
	friend
	inline
	xyzVector
	center( xyzVector const & a, xyzVector const & b, xyzVector const & c, xyzVector const & d )
	{
		return xyzVector(
			Value( 0.25 * ( a.x_ + b.x_ + c.x_ + d.x_ ) ),
			Value( 0.25 * ( a.y_ + b.y_ + c.y_ + d.y_ ) ),
			Value( 0.25 * ( a.z_ + b.z_ + c.z_ + d.z_ ) )
		);
	}


	/// @brief Center of 4 xyzVectors: Return via argument (slightly faster)
	friend
	inline
	void
	center( xyzVector const & a, xyzVector const & b, xyzVector const & c, xyzVector const & d, xyzVector & m )
	{
		m.x_ = Value( 0.25 * ( a.x_ + b.x_ + c.x_ + d.x_ ) );
		m.y_ = Value( 0.25 * ( a.y_ + b.y_ + c.y_ + d.y_ ) );
		m.z_ = Value( 0.25 * ( a.z_ + b.z_ + c.z_ + d.z_ ) );
	}


	/// @brief Angle between two vectors (in radians on [ 0, pi ])
	friend
	inline
	Value
	angle_of( xyzVector const & a, xyzVector const & b )
	{
		Value const mag = a.length() * b.length();
		return ( mag > Value( 0 ) ? std::acos( sin_cos_range( a.dot( b ) / mag ) ) : Value( 0 ) );
	}


	/// @brief Angle formed by three consecutive points (in radians on [ 0, pi ])
	/// @note  For points a, b, c, the angle is the angle between the vectors a - b  and c - b
	///        in other words, the positive angle about b from a to c
	friend
	inline
	Value
	angle_of( xyzVector const & a, xyzVector const & b, xyzVector const & c )
	{
		return angle_of( a - b, c - b );
	}


	/// @brief Cosine of angle between two vectors
	friend
	inline
	Value
	cos_of( xyzVector const & a, xyzVector const & b )
	{
		Value const mag = a.length() * b.length();
		return ( mag > Value( 0 ) ? sin_cos_range( a.dot( b ) / mag ) : Value( 1 ) );
	}


	/// @brief Cosine of angle formed by three consecutive points
	/// @note  For points a, b, c, the angle is the angle between the vectors a - b  and c - b
	///        in other words, the positive angle about b from a to c.
	friend
	inline
	Value
	cos_of( xyzVector const & a, xyzVector const & b, xyzVector const & c )
	{
		return cos_of( a - b, c - b );
	}


	/// @brief Sine of angle between two vectors
	friend
	inline
	Value
	sin_of( xyzVector const & a, xyzVector const & b )
	{
		return std::sqrt( Value( 1 ) - square( cos_of( a, b ) ) );
	}


	/// @brief Sine of angle formed by three consecutive points
	/// @note  For points a, b, c, the angle is the angle between the vectors a - b  and c - b
	///        in other words, the positive angle about b from a to c
	friend
	inline
	Value
	sin_of( xyzVector const & a, xyzVector const & b, xyzVector const & c )
	{
		return sin_of( a - b, c - b );
	}

	// AMW: for md code
	// TODO: figure out what this is in fundamental vector operations
	friend
	inline
	xyzVector
	update_operation( xyzVector const & a, xyzVector const & b )
	{
		return xyzVector(
						   a.y_ * a.y_ * b.x_
						 + a.z_ * a.z_ * b.x_
						 - a.x_ * a.y_ * b.y_
						 - a.x_ * a.z_ * b.z_,
						 
						   a.z_ * a.z_ * b.y_
						 + a.x_ * a.x_ * b.y_
						 - a.y_ * a.z_ * b.z_
						 - a.y_ * a.x_ * b.x_,
						 
						   a.x_ * a.x_ * b.z_
						 + a.y_ * a.y_ * b.z_
						 - a.z_ * a.x_ * b.x_
						 - a.z_ * a.y_ * b.y_
						 
						 );
	}
	
	// AMW: for md code
	// TODO: figure out what this is in fundamental vector operations
	friend
	inline
	xyzVector
	update_5way_operation(
		xyzVector const & a,
		xyzVector const & b,
		xyzVector const & c,
		xyzVector const & d,
		xyzVector const & e
	) {
		return xyzVector(
						 - a.y_ * b.y_ * c.x_
						 - a.z_ * b.z_ * c.x_
						 + a.x_ * b.y_ * c.y_
						 + a.x_ * b.y_ * c.y_
						 - b.x_ * a.y_ * c.y_
						 + a.x_ * b.z_ * c.z_
						 + a.x_ * b.z_ * c.z_
						 - b.x_ * a.z_ * c.z_
						 + d.z_ * e.y_
						 - d.y_ * e.z_,
						 
						 - a.z_ * b.z_ * c.y_
						 - a.x_ * b.x_ * c.y_
						 + a.y_ * b.z_ * c.z_
						 + a.y_ * b.z_ * c.z_
						 - b.y_ * a.z_ * c.z_
						 + a.y_ * b.x_ * c.x_
						 + a.y_ * b.x_ * c.x_
						 - b.y_ * a.x_ * c.x_
						 + d.x_ * e.z_
						 - d.z_ * e.x_,
						 
						 - a.x_ * b.x_ * c.z_
						 - a.y_ * b.y_ * c.z_
						 + a.z_ * b.x_ * c.x_
						 + a.z_ * b.x_ * c.x_
						 - b.z_ * a.x_ * c.x_
						 + a.z_ * b.x_ * c.y_
						 + a.z_ * b.x_ * c.y_
						 - b.z_ * a.y_ * c.y_
						 + d.y_ * e.x_
						 - d.x_ * e.y_
						 
						 );
	}

public: // Properties: predicates


	/// @brief Is zero?
	inline
	bool
	is_zero() const
	{
		static Value const ZERO( 0 );
		return ( x_ == ZERO ) && ( y_ == ZERO ) && ( z_ == ZERO );
	}


	/// @brief Is exactly normalized?
	inline
	bool
	is_normalized() const
	{
		return ( length_squared() == Value( 1 ) );
	}


	/// @brief Is normalized to within a tolerance?
	inline
	bool
	is_normalized( Value const & tol ) const
	{
		Value const tol_sq = tol * tol;
		Value const length_sq = length_squared();
		return ( ( length_sq >= Value( 1 ) - tol_sq ) && ( length_sq <= Value( 1 ) + tol_sq ) );
	}


	/// @brief Is exactly a unit vector?
	inline
	bool
	is_unit() const
	{
		return ( length_squared() == Value( 1 ) );
	}


	/// @brief Is a unit vector to within a tolerance?
	inline
	bool
	is_unit( Value const & tol ) const
	{
		Value const tol_sq = tol * tol;
		Value const length_sq = length_squared();
		return ( ( length_sq >= Value( 1 ) - tol_sq ) && ( length_sq <= Value( 1 ) + tol_sq ) );
	}


public: // Properties: accessors


	/// @brief Value x const
	inline
	Value const &
	x() const
	{
		return x_;
	}


	/// @brief Value x
	inline
	Value &
	x()
	{
		return x_;
	}


	/// @brief Value y const
	inline
	Value const &
	y() const
	{
		return y_;
	}


	/// @brief Value y
	inline
	Value &
	y()
	{
		return y_;
	}


	/// @brief Value z const
	inline
	Value const &
	z() const
	{
		return z_;
	}


	/// @brief Value z
	inline
	Value &
	z()
	{
		return z_;
	}


	/// @brief Length
	inline
	Value
	length() const
	{
		return std::sqrt( ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Length squared
	inline
	Value
	length_squared() const
	{
		return ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ );
	}


	/// @brief Norm
	inline
	Value
	norm() const
	{
		return std::sqrt( ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Norm squared
	inline
	Value
	norm_squared() const
	{
		return ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ );
	}


	/// @brief Magnitude
	inline
	Value
	magnitude() const
	{
		return std::sqrt( ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ ) );
	}


	/// @brief Magnitude squared
	inline
	Value
	magnitude_squared() const
	{
		return ( x_ * x_ ) + ( y_ * y_ ) + ( z_ * z_ );
	}


public: // Properties: value assignment


	/// @brief x assignment
	inline
	void
	x( Value const & x_a )
	{
		x_ = x_a;
	}


	/// @brief y assignment
	inline
	void
	y( Value const & y_a )
	{
		y_ = y_a;
	}


	/// @brief z assignment
	inline
	void
	z( Value const & z_a )
	{
		z_ = z_a;
	}


public: // Indexers

	/// @brief xyzVector.at: 0-based index with bounds checking
	inline
	Value const &
	at(int const i) const
	{
		if(!((i >= 0) && (i < 3)))
		{
			throw std::out_of_range ("numeric::xyzVector::at");
		}

		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}

	/// @brief xyzVector.at: 0-based index with bounds checking
	inline
	Value &
	at(int const i)
	{
		if(!((i >= 0) && (i < 3)))
		{
			throw std::out_of_range ("numeric::xyzVector::at");
		}

		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}

	/// @brief xyzVector[ i ] const: 0-based index
	inline
	Value const &
	operator []( int const i ) const
	{
		assert( ( i >= 0 ) && ( i < 3 ) );
		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}


	/// @brief xyzVector[ i ]: 0-based index
	inline
	Value &
	operator []( int const i )
	{
		assert( ( i >= 0 ) && ( i < 3 ) );
		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}


	/// @brief xyzVector( i ) const: 1-based index
	inline
	Value const &
	operator ()( int const i ) const
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		return ( i == 1 ? x_ : ( i == 2 ? y_ : z_ ) );
	}


	/// @brief xyzVector( i ): 1-based index
	inline
	Value &
	operator ()( int const i )
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		return ( i == 1 ? x_ : ( i == 2 ? y_ : z_ ) );
	}


public: // Comparison


	/// @brief xyzVector == xyzVector
	friend
	inline
	bool
	operator ==( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ == b.x_ ) && ( a.y_ == b.y_ ) && ( a.z_ == b.z_ );
	}


	/// @brief xyzVector != xyzVector
	friend
	inline
	bool
	operator !=( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ != b.x_ ) || ( a.y_ != b.y_ ) || ( a.z_ != b.z_ );
	}


	/// @brief xyzVector < xyzVector
	friend
	inline
	bool
	operator <( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ < b.x_ ) && ( a.y_ < b.y_ ) && ( a.z_ < b.z_ );
	}


	/// @brief xyzVector <= xyzVector
	friend
	inline
	bool
	operator <=( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ <= b.x_ ) && ( a.y_ <= b.y_ ) && ( a.z_ <= b.z_ );
	}


	/// @brief xyzVector >= xyzVector
	friend
	inline
	bool
	operator >=( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ >= b.x_ ) && ( a.y_ >= b.y_ ) && ( a.z_ >= b.z_ );
	}


	/// @brief xyzVector > xyzVector
	friend
	inline
	bool
	operator >( xyzVector const & a, xyzVector const & b )
	{
		return ( a.x_ > b.x_ ) && ( a.y_ > b.y_ ) && ( a.z_ > b.z_ );
	}


	/// @brief xyzVector == Value
	friend
	inline
	bool
	operator ==( xyzVector const & v, Value const & t )
	{
		return ( v.x_ == t ) && ( v.y_ == t ) && ( v.z_ == t );
	}


	/// @brief xyzVector != Value
	friend
	inline
	bool
	operator !=( xyzVector const & v, Value const & t )
	{
		return ( v.x_ != t ) || ( v.y_ != t ) || ( v.z_ != t );
	}


	/// @brief xyzVector < Value
	friend
	inline
	bool
	operator <( xyzVector const & v, Value const & t )
	{
		return ( v.x_ < t ) && ( v.y_ < t ) && ( v.z_ < t );
	}


	/// @brief xyzVector <= Value
	friend
	inline
	bool
	operator <=( xyzVector const & v, Value const & t )
	{
		return ( v.x_ <= t ) && ( v.y_ <= t ) && ( v.z_ <= t );
	}


	/// @brief xyzVector >= Value
	friend
	inline
	bool
	operator >=( xyzVector const & v, Value const & t )
	{
		return ( v.x_ >= t ) && ( v.y_ >= t ) && ( v.z_ >= t );
	}


	/// @brief xyzVector > Value
	friend
	inline
	bool
	operator >( xyzVector const & v, Value const & t )
	{
		return ( v.x_ > t ) && ( v.y_ > t ) && ( v.z_ > t );
	}


	/// @brief Value == xyzVector
	friend
	inline
	bool
	operator ==( Value const & t, xyzVector const & v )
	{
		return ( t == v.x_ ) && ( t == v.y_ ) && ( t == v.z_ );
	}


	/// @brief Value != xyzVector
	friend
	inline
	bool
	operator !=( Value const & t, xyzVector const & v )
	{
		return ( t != v.x_ ) || ( t != v.y_ ) || ( t != v.z_ );
	}


	/// @brief Value < xyzVector
	friend
	inline
	bool
	operator <( Value const & t, xyzVector const & v )
	{
		return ( t < v.x_ ) && ( t < v.y_ ) && ( t < v.z_ );
	}


	/// @brief Value <= xyzVector
	friend
	inline
	bool
	operator <=( Value const & t, xyzVector const & v )
	{
		return ( t <= v.x_ ) && ( t <= v.y_ ) && ( t <= v.z_ );
	}


	/// @brief Value >= xyzVector
	friend
	inline
	bool
	operator >=( Value const & t, xyzVector const & v )
	{
		return ( t >= v.x_ ) && ( t >= v.y_ ) && ( t >= v.z_ );
	}


	/// @brief Value > xyzVector
	friend
	inline
	bool
	operator >( Value const & t, xyzVector const & v )
	{
		return ( t > v.x_ ) && ( t > v.y_ ) && ( t > v.z_ );
	}


	/// @brief Equal length?
	inline
	bool
	equal_length( xyzVector const & v )
	{
		return ( length_squared() == v.length_squared() );
	}


	/// @brief Equal length?
	friend
	inline
	bool
	equal_length( xyzVector const & a, xyzVector const & b )
	{
		return ( a.length_squared() == b.length_squared() );
	}


	/// @brief Not equal length?
	inline
	bool
	not_equal_length( xyzVector const & v )
	{
		return ( length_squared() != v.length_squared() );
	}


	/// @brief Not equal length?
	friend
	inline
	bool
	not_equal_length( xyzVector const & a, xyzVector const & b )
	{
		return ( a.length_squared() != b.length_squared() );
	}


	/// @brief Longer?
	inline
	bool
	longer( xyzVector const & v )
	{
		return ( length_squared() > v.length_squared() );
	}


	/// @brief Longer or equal length?
	inline
	bool
	longer_or_equal( xyzVector const & v )
	{
		return ( length_squared() >= v.length_squared() );
	}


	/// @brief Shorter?
	inline
	bool
	shorter( xyzVector const & v )
	{
		return ( length_squared() < v.length_squared() );
	}


	/// @brief Shorter or equal length?
	inline
	bool
	shorter_or_equal( xyzVector const & v )
	{
		return ( length_squared() <= v.length_squared() );
	}

	/// @brief Hashing of coords using boost::hash
	template <typename U>
	friend platform::Size hash_value(xyzVector<U> const & v);

private: // Methods


	/// @brief square( t ) == t * t
	inline
	static
	Value
	square( Value const & t )
	{
		return t * t;
	}


private: // Fields


	/// @brief Coordinates of the 3 coordinate vector
	Value x_;
	Value y_;
	Value z_;


}; // xyzVector


//// @brief Hashing of coords using boost::hash
template< typename T >
platform::Size hash_value(xyzVector< T >  const & v)
{
	std::size_t hash = 0;
	boost::hash_combine(hash, v.x_);
	boost::hash_combine(hash, v.y_);
	boost::hash_combine(hash, v.z_);
	return hash;
}


/// @brief xyzVector + xyzVector
template< typename T >
xyzVector< T >
operator +( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector + T
template< typename T >
xyzVector< T >
operator +( xyzVector< T > const & v, T const & t );


/// @brief T + xyzVector
template< typename T >
xyzVector< T >
operator +( T const & t, xyzVector< T > const & v );


/// @brief xyzVector - xyzVector
template< typename T >
xyzVector< T >
operator -( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector - T
template< typename T >
xyzVector< T >
operator -( xyzVector< T > const & v, T const & t );


/// @brief T - xyzVector
template< typename T >
xyzVector< T >
operator -( T const & t, xyzVector< T > const & v );


/// @brief xyzVector * T
template< typename T >
xyzVector< T >
operator *( xyzVector< T > const & v, T const & t );


/// @brief T * xyzVector
template< typename T >
xyzVector< T >
operator *( T const & t, xyzVector< T > const & v );


/// @brief xyzVector / T
template< typename T >
xyzVector< T >
operator /( xyzVector< T > const & v, T const & t );


/// @brief Add: xyzVector + xyzVector
template< typename T >
void
add( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > & r );


/// @brief Add: xyzVector + T
template< typename T >
void
add( xyzVector< T > const & v, T const & t, xyzVector< T > & r );


/// @brief Add: T + xyzVector
template< typename T >
void
add( T const & t, xyzVector< T > const & v, xyzVector< T > & r );


/// @brief Subtract: xyzVector - xyzVector
template< typename T >
void
subtract( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > & r );


/// @brief Subtract: xyzVector - T
template< typename T >
void
subtract( xyzVector< T > const & v, T const & t, xyzVector< T > & r );


/// @brief Subtract: T - xyzVector
template< typename T >
void
subtract( T const & t, xyzVector< T > const & v, xyzVector< T > & r );


/// @brief Multiply: xyzVector * T
template< typename T >
void
multiply( xyzVector< T > const & v, T const & t, xyzVector< T > & r );


/// @brief Multiply: T * xyzVector
template< typename T >
void
multiply( T const & t, xyzVector< T > const & v, xyzVector< T > & r );


/// @brief Divide: xyzVector / T
template< typename T >
void
divide( xyzVector< T > const & v, T const & t, xyzVector< T > & r );


/// @brief xyzVector with min coordinates of two xyzVectors
template< typename T >
xyzVector< T >
min( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector with max coordinates of two xyzVectors
template< typename T >
xyzVector< T >
max( xyzVector< T > const & a, xyzVector< T > const & b );

//Commented out 8/9/11 SML; xyzVector friend functions cause problems for some versions of GCC.  Read the mailing list logs for 2/17/09 for further information.
/// @brief Distance
// template< typename T >
// T
// distance( xyzVector< T > const & a, xyzVector< T > const & b );

//Commented out 8/9/11 SML; xyzVector friend functions cause problems for some versions of GCC.  Read the mailing list logs for 2/17/09 for further information.
/// @brief Distance squared
// template< typename T >
// T
// distance_squared( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Dot product
template< typename T >
T
dot( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Dot product
template< typename T >
T
dot_product( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Inner product ( == dot product )
template< typename T >
T
inner_product( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Cross product
template< typename T >
xyzVector< T >
cross( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Cross product
template< typename T >
xyzVector< T >
cross_product( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Cross product: Return via argument (slightly faster)
template< typename T >
void
cross( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > & c );


/// @brief Cross product: Return via argument (slightly faster)
template< typename T >
void
cross_product( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > & c );


/// @brief Midpoint of 2 xyzVectors
template< typename T >
xyzVector< T >
midpoint( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Midpoint of 2 xyzVectors: Return via argument (slightly faster)
template< typename T >
void
midpoint( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > & m );


/// @brief Center of 2 xyzVectors
template< typename T >
xyzVector< T >
center( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Center of 2 xyzVectors: Return via argument (slightly faster)
template< typename T >
void
center( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > & m );


/// @brief Center of 3 xyzVectors
template< typename T >
xyzVector< T >
center( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c );


/// @brief Center of 3 xyzVectors: Return via argument (slightly faster)
template< typename T >
void
center( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c, xyzVector< T > & m );


/// @brief Center of 4 xyzVectors
template< typename T >
xyzVector< T >
center( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c, xyzVector< T > const & d );


/// @brief Center of 4 xyzVectors: Return via argument (slightly faster)
template< typename T >
void
center( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c, xyzVector< T > const & d, xyzVector< T > & m );


/// @brief Angle between two vectors (in radians on [ 0, pi ])
template< typename T >
T
angle_of( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Angle formed by three consecutive points (in radians on [ 0, pi ])
template< typename T >
T
angle_of( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c );


/// @brief Cosine of angle between two vectors
template< typename T >
T
cos_of( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Cosine of angle formed by three consecutive points
template< typename T >
T
cos_of( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c );


/// @brief Sine of angle between two vectors
template< typename T >
T
sin_of( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Sine of angle formed by three consecutive points
template< typename T >
T
sin_of( xyzVector< T > const & a, xyzVector< T > const & b, xyzVector< T > const & c );


/// @brief xyzVector == xyzVector
template< typename T >
bool
operator ==( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector != xyzVector
template< typename T >
bool
operator !=( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector < xyzVector
template< typename T >
bool
operator <( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector <= xyzVector
template< typename T >
bool
operator <=( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector >= xyzVector
template< typename T >
bool
operator >=( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector > xyzVector
template< typename T >
bool
operator >( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief xyzVector == T
template< typename T >
bool
operator ==( xyzVector< T > const & v, T const & t );


/// @brief xyzVector != T
template< typename T >
bool
operator !=( xyzVector< T > const & v, T const & t );


/// @brief xyzVector < T
template< typename T >
bool
operator <( xyzVector< T > const & v, T const & t );


/// @brief xyzVector <= T
template< typename T >
bool
operator <=( xyzVector< T > const & v, T const & t );


/// @brief xyzVector >= T
template< typename T >
bool
operator >=( xyzVector< T > const & v, T const & t );


/// @brief xyzVector > T
template< typename T >
bool
operator >( xyzVector< T > const & v, T const & t );


/// @brief T == xyzVector
template< typename T >
bool
operator ==( T const & t, xyzVector< T > const & v );


/// @brief T != xyzVector
template< typename T >
bool
operator !=( T const & t, xyzVector< T > const & v );


/// @brief T < xyzVector
template< typename T >
bool
operator <( T const & t, xyzVector< T > const & v );


/// @brief T <= xyzVector
template< typename T >
bool
operator <=( T const & t, xyzVector< T > const & v );


/// @brief T >= xyzVector
template< typename T >
bool
operator >=( T const & t, xyzVector< T > const & v );


/// @brief T > xyzVector
template< typename T >
bool
operator >( T const & t, xyzVector< T > const & v );


/// @brief Equal length?
template< typename T >
bool
equal_length( xyzVector< T > const & a, xyzVector< T > const & b );


/// @brief Not equal length?
template< typename T >
bool
not_equal_length( xyzVector< T > const & a, xyzVector< T > const & b );


// PyRosetta work around for templates classes
//class xyzVector_Double : public xyzVector< double >
//{};

// we don't need to define it here, it already defined in numeric/xyzVector.io.hh
// template< typename T >
// std::ostream & operator <<(std::ostream & os, xyzVector< T > const & v) { os << "[x=" << v.x() << ", y=" << v.y() << ", z=" << v.z() << "]"; }


} // namespace numeric

namespace ObjexxFCL {

/// @brief Specialization for FArrayTraits, to allow reasonable default constructor in FArray context
template< typename T >
struct FArrayTraits < numeric::xyzVector< T > >
{
	typedef numeric::xyzVector< T >  traits_type;

	/// @brief Initial Value
	inline
	static
	traits_type
	initial_value()
	{
	        return numeric::xyzVector< T >(0,0,0); // Use all zeros
	}
}; // FArrayTraits

// stream operator for show function
template< typename T >
std::ostream &operator<< (std::ostream &os, numeric::xyzVector< T > const &vector)
{
	vector.show(os);
	return os;
}

}


#endif // INCLUDED_numeric_xyzVector_HH
