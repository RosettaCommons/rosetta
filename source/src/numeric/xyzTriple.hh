// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/xyzTriple.hh
/// @brief  Fast (x,y,z)-coordinate vector container
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Inline, loop-free functions for speed
///  @li Non-virtual destructor for speed: Not set up for use as a base class
///  @li Pointer constructor and assignment not available for xyzTriples of pointers
///  @li Container semantics: lexicographic complete ordering


#ifndef INCLUDED_numeric_xyzTriple_hh
#define INCLUDED_numeric_xyzTriple_hh


// Unit headers
#include <numeric/xyzTriple.fwd.hh>

// Package headers
#include <numeric/trig.functions.hh>

// C++ headers
#include <utility/assert.hh>
#include <cmath>
#include <stdexcept>


namespace numeric {


/// @brief Fast (x,y,z)-coordinate vector container
template< typename T >
class xyzTriple
{


private: // Friends


	template< typename > friend class xyzTriple;


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
	xyzTriple()
	{}


	/// @brief Copy constructor
	inline
	xyzTriple( xyzTriple const & v ) :
		x_( v.x_ ),
		y_( v.y_ ),
		z_( v.z_ )
	{}


	/// @brief Copy constructor
	template< typename U >
	inline
	xyzTriple( xyzTriple< U > const & v ) :
		x_( v.x_ ),
		y_( v.y_ ),
		z_( v.z_ )
	{}


	/// @brief Uniform value constructor
	inline
	explicit
	xyzTriple( Value const & t ) :
		x_( t ),
		y_( t ),
		z_( t )
	{}


	/// @brief Triple value constructor
	inline
	xyzTriple(
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
	xyzTriple( U const * p ) :
		x_( *p ),
		y_( *++p ),
		z_( *++p )
	{}


	/// @brief Destructor
	inline
	~xyzTriple()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	xyzTriple &
	operator =( xyzTriple const & v )
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
	xyzTriple &
	operator =( xyzTriple< U > const & v )
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
	xyzTriple &
	operator =( U const * p )
	{
		x_ = *p;
		y_ = *++p;
		z_ = *++p;
		return *this;
	}


	/// @brief += xyzTriple
	template< typename U >
	inline
	xyzTriple &
	operator +=( xyzTriple< U > const & v )
	{
		x_ += v.x_;
		y_ += v.y_;
		z_ += v.z_;
		return *this;
	}


	/// @brief -= xyzTriple
	template< typename U >
	inline
	xyzTriple &
	operator -=( xyzTriple< U > const & v )
	{
		x_ -= v.x_;
		y_ -= v.y_;
		z_ -= v.z_;
		return *this;
	}


	/// @brief Assign Value * xyzTriple
	/// @note  Avoids temporary of = Value * xyzTriple
	template< typename U >
	inline
	xyzTriple &
	scaled_assign( Value const & t, xyzTriple< U > const & v )
	{
		x_ = t * v.x_;
		y_ = t * v.y_;
		z_ = t * v.z_;
		return *this;
	}


	/// @brief Add Value * xyzTriple
	/// @note  Avoids temporary of += Value * xyzTriple
	template< typename U >
	inline
	xyzTriple &
	scaled_add( Value const & t, xyzTriple< U > const & v )
	{
		x_ += t * v.x_;
		y_ += t * v.y_;
		z_ += t * v.z_;
		return *this;
	}


	/// @brief Subtract Value * xyzTriple
	/// @note  Avoids temporary of -= Value * xyzTriple
	template< typename U >
	inline
	xyzTriple &
	scaled_sub( Value const & t, xyzTriple< U > const & v )
	{
		x_ -= t * v.x_;
		y_ -= t * v.y_;
		z_ -= t * v.z_;
		return *this;
	}


	/// @brief = Value
	inline
	xyzTriple &
	operator =( Value const & t )
	{
		x_ = y_ = z_ = t;
		return *this;
	}


	/// @brief += Value
	inline
	xyzTriple &
	operator +=( Value const & t )
	{
		x_ += t;
		y_ += t;
		z_ += t;
		return *this;
	}


	/// @brief -= Value
	inline
	xyzTriple &
	operator -=( Value const & t )
	{
		x_ -= t;
		y_ -= t;
		z_ -= t;
		return *this;
	}


	/// @brief *= Value
	inline
	xyzTriple &
	operator *=( Value const & t )
	{
		x_ *= t;
		y_ *= t;
		z_ *= t;
		return *this;
	}


	/// @brief /= Value
	inline
	xyzTriple &
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
	xyzTriple &
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


	/// @brief Clear
	inline
	xyzTriple &
	clear()
	{
		x_ = y_ = z_ = Value( 0 );
		return *this;
	}


	/// @brief Zero
	inline
	xyzTriple &
	zero()
	{
		x_ = y_ = z_ = Value( 0 );
		return *this;
	}


	/// @brief Negate
	inline
	xyzTriple &
	negate()
	{
		x_ = -x_;
		y_ = -y_;
		z_ = -z_;
		return *this;
	}


	/// @brief -xyzTriple (negated copy)
	inline
	xyzTriple
	operator -() const
	{
		return xyzTriple( -x_, -y_, -z_ );
	}


	/// @brief Negated copy
	inline
	xyzTriple
	negated() const
	{
		return xyzTriple( -x_, -y_, -z_ );
	}


	/// @brief Negated: Return via argument (slightly faster)
	inline
	void
	negated( xyzTriple & a ) const
	{
		a.x_ = -x_;
		a.y_ = -y_;
		a.z_ = -z_;
	}


	/// @brief xyzTriple + xyzTriple
	friend
	inline
	xyzTriple
	operator +( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple( a.x_ + b.x_, a.y_ + b.y_, a.z_ + b.z_ );
	}


	/// @brief xyzTriple + Value
	friend
	inline
	xyzTriple
	operator +( xyzTriple const & v, Value const & t )
	{
		return xyzTriple( v.x_ + t, v.y_ + t, v.z_ + t );
	}


	/// @brief Value + xyzTriple
	friend
	inline
	xyzTriple
	operator +( Value const & t, xyzTriple const & v )
	{
		return xyzTriple( t + v.x_, t + v.y_, t + v.z_ );
	}


	/// @brief xyzTriple - xyzTriple
	friend
	inline
	xyzTriple
	operator -( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple( a.x_ - b.x_, a.y_ - b.y_, a.z_ - b.z_ );
	}


	/// @brief xyzTriple - Value
	friend
	inline
	xyzTriple
	operator -( xyzTriple const & v, Value const & t )
	{
		return xyzTriple( v.x_ - t, v.y_ - t, v.z_ - t );
	}


	/// @brief Value - xyzTriple
	friend
	inline
	xyzTriple
	operator -( Value const & t, xyzTriple const & v )
	{
		return xyzTriple( t - v.x_, t - v.y_, t - v.z_ );
	}


	/// @brief xyzTriple * Value
	friend
	inline
	xyzTriple
	operator *( xyzTriple const & v, Value const & t )
	{
		return xyzTriple( v.x_ * t, v.y_ * t, v.z_ * t );
	}


	/// @brief Value * xyzTriple
	friend
	inline
	xyzTriple
	operator *( Value const & t, xyzTriple const & v )
	{
		return xyzTriple( t * v.x_, t * v.y_, t * v.z_ );
	}


	/// @brief xyzTriple / Value
	friend
	inline
	xyzTriple
	operator /( xyzTriple const & v, Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value ( 1 ) / t );
		return xyzTriple( v.x_ * inv_t, v.y_ * inv_t, v.z_ * inv_t );
	}


	/// @brief Add: xyzTriple + xyzTriple
	friend
	inline
	void
	add( xyzTriple const & a, xyzTriple const & b, xyzTriple & r )
	{
		r.x_ = a.x_ + b.x_;
		r.y_ = a.y_ + b.y_;
		r.z_ = a.z_ + b.z_;
	}


	/// @brief Add: xyzTriple + Value
	friend
	inline
	void
	add( xyzTriple const & v, Value const & t, xyzTriple & r )
	{
		r.x_ = v.x_ + t;
		r.y_ = v.y_ + t;
		r.z_ = v.z_ + t;
	}


	/// @brief Add: Value + xyzTriple
	friend
	inline
	void
	add( Value const & t, xyzTriple const & v, xyzTriple & r )
	{
		r.x_ = t + v.x_;
		r.y_ = t + v.y_;
		r.z_ = t + v.z_;
	}


	/// @brief Subtract: xyzTriple - xyzTriple
	friend
	inline
	void
	subtract( xyzTriple const & a, xyzTriple const & b, xyzTriple & r )
	{
		r.x_ = a.x_ - b.x_;
		r.y_ = a.y_ - b.y_;
		r.z_ = a.z_ - b.z_;
	}


	/// @brief Subtract: xyzTriple - Value
	friend
	inline
	void
	subtract( xyzTriple const & v, Value const & t, xyzTriple & r )
	{
		r.x_ = v.x_ - t;
		r.y_ = v.y_ - t;
		r.z_ = v.z_ - t;
	}


	/// @brief Subtract: Value - xyzTriple
	friend
	inline
	void
	subtract( Value const & t, xyzTriple const & v, xyzTriple & r )
	{
		r.x_ = t - v.x_;
		r.y_ = t - v.y_;
		r.z_ = t - v.z_;
	}


	/// @brief Multiply: xyzTriple * Value
	friend
	inline
	void
	multiply( xyzTriple const & v, Value const & t, xyzTriple & r )
	{
		r.x_ = v.x_ * t;
		r.y_ = v.y_ * t;
		r.z_ = v.z_ * t;
	}


	/// @brief Multiply: Value * xyzTriple
	friend
	inline
	void
	multiply( Value const & t, xyzTriple const & v, xyzTriple & r )
	{
		r.x_ = t * v.x_;
		r.y_ = t * v.y_;
		r.z_ = t * v.z_;
	}


	/// @brief Divide: xyzTriple / Value
	friend
	inline
	void
	divide( xyzTriple const & v, Value const & t, xyzTriple & r )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		r.x_ = v.x_ * inv_t;
		r.y_ = v.y_ * inv_t;
		r.z_ = v.z_ * inv_t;
	}


	/// @brief Set minimum coordinates wrt another xyzTriple
	inline
	xyzTriple &
	min( xyzTriple const & v )
	{
		x_ = ( x_ <= v.x_ ? x_ : v.x_ );
		y_ = ( y_ <= v.y_ ? y_ : v.y_ );
		z_ = ( z_ <= v.z_ ? z_ : v.z_ );
		return *this;
	}


	/// @brief Set maximum coordinates wrt another xyzTriple
	inline
	xyzTriple &
	max( xyzTriple const & v )
	{
		x_ = ( x_ >= v.x_ ? x_ : v.x_ );
		y_ = ( y_ >= v.y_ ? y_ : v.y_ );
		z_ = ( z_ >= v.z_ ? z_ : v.z_ );
		return *this;
	}


	/// @brief xyzTriple with min coordinates of two xyzTriples
	friend
	inline
	xyzTriple
	min( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple(
		 ( a.x_ <= b.x_ ? a.x_ : b.x_ ),
		 ( a.y_ <= b.y_ ? a.y_ : b.y_ ),
		 ( a.z_ <= b.z_ ? a.z_ : b.z_ )
		);
	}


	/// @brief xyzTriple with max coordinates of two xyzTriples
	friend
	inline
	xyzTriple
	max( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple(
		 ( a.x_ >= b.x_ ? a.x_ : b.x_ ),
		 ( a.y_ >= b.y_ ? a.y_ : b.y_ ),
		 ( a.z_ >= b.z_ ? a.z_ : b.z_ )
		);
	}


	/// @brief Normalize
	inline
	xyzTriple &
	normalize()
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const inv_length( Value( 1 ) / length_ );
		x_ *= inv_length;
		y_ *= inv_length;
		z_ *= inv_length;
		return *this;
	}


	/// @brief Normalized
	inline
	void
	normalized( xyzTriple & a ) const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const inv_length( Value( 1 ) / length_ );
		a.x_ = x_ * inv_length;
		a.y_ = y_ * inv_length;
		a.z_ = z_ * inv_length;
	}


	/// @brief Normalize: zero xyzTriple if length is zero
	inline
	xyzTriple &
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


	/// @brief Normalized: zero xyzTriple if length is zero
	inline
	void
	normalized_or_zero( xyzTriple & a ) const
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


	/// @brief Normalize: arbitrary normalized xyzTriple if length is zero
	inline
	xyzTriple &
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


	/// @brief Normalized: arbitrary normalized xyzTriple if length is zero
	inline
	void
	normalized_any( xyzTriple & a ) const
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
	xyzTriple &
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
	normalized( Value const & length_a, xyzTriple & a ) const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const dilation = length_a / length_;
		a.x_ = x_ * dilation;
		a.y_ = y_ * dilation;
		a.z_ = z_ * dilation;
	}


	/// @brief Normalize to a length: zero xyzTriple if length is zero
	inline
	xyzTriple &
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


	/// @brief Normalized to a length: zero xyzTriple if length is zero
	inline
	void
	normalized_or_zero( Value const & length_a, xyzTriple & a ) const
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


	/// @brief Normalize to a length: arbitrary normalized xyzTriple if length is zero
	inline
	xyzTriple &
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


	/// @brief Normalized to a length: arbitrary normalized xyzTriple if length is zero
	inline
	void
	normalized_any( Value const & length_a, xyzTriple & a ) const
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
	xyzTriple
	normalized() const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const inv_length( Value( 1 ) / length_ );
		return xyzTriple( x_ * inv_length, y_ * inv_length, z_ * inv_length );
	}


	/// @brief Normalized copy: Zero xyzTriple if length is zero
	inline
	xyzTriple
	normalized_or_zero() const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			return xyzTriple( x_ * inv_length, y_ * inv_length, z_ * inv_length );
		} else { // Return zero vector
			return xyzTriple( Value( 0 ), Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Normalized copy: Arbitrary normalized xyzTriple if length is zero
	inline
	xyzTriple
	normalized_any() const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const inv_length( Value( 1 ) / length_ );
			return xyzTriple( x_ * inv_length, y_ * inv_length, z_ * inv_length );
		} else { // Return arbitrary normalized vector
			return xyzTriple( Value( 1 ), Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Normalized to a length copy
	inline
	xyzTriple
	normalized( Value const & length_a ) const
	{
		Value const length_ = length();
		assert( length_ != Value ( 0 ) );
		Value const dilation = length_a / length_;
		return xyzTriple( x_ * dilation, y_ * dilation, z_ * dilation );
	}


	/// @brief Normalized to a length copy: Zero xyzTriple if length is zero
	inline
	xyzTriple
	normalized_or_zero( Value const & length_a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			return xyzTriple( x_ * dilation, y_ * dilation, z_ * dilation );
		} else { // Return zero vector
			return xyzTriple( Value( 0 ), Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Normalized to a length copy: Arbitrary normalized xyzTriple if length is zero
	inline
	xyzTriple
	normalized_any( Value const & length_a ) const
	{
		Value const length_ = length();
		if ( length_ > Value( 0 ) ) {
			Value const dilation = length_a / length_;
			return xyzTriple( x_ * dilation, y_ * dilation, z_ * dilation );
		} else { // Return arbitrary normalized vector
			return xyzTriple( length_a, Value( 0 ), Value( 0 ) );
		}
	}


	/// @brief Project normal
	/// @note  This vector projected normally onto input vector
	/// @note  Not meaningful when v == 0
	inline
	xyzTriple &
	project_normal( xyzTriple const & v )
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
	xyzTriple
	projected_normal( xyzTriple const & v ) const
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		return xyzTriple( x_ - ( c * v.x_ ), y_ - ( c * v.y_ ), z_ - ( c * v.z_ ) );
	}


	/// @brief Projected normal
	/// @note  Copy of this vector projected normally onto first input vector
	/// @note  Not meaningful when v == 0
	inline
	void
	projected_normal( xyzTriple const & v, xyzTriple & a ) const
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
	xyzTriple &
	project_parallel( xyzTriple const & v )
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
	xyzTriple
	projected_parallel( xyzTriple const & v ) const
	{
		assert( v.length_squared() != Value( 0 ) );
		Value const c = dot( v ) / v.length_squared();
		return xyzTriple( c * v.x_, c * v.y_, c * v.z_ );
	}


	/// @brief Projected parallel
	/// @note  Copy of this vector projected in direction of first input vector
	/// @note  Not meaningful when v == 0
	inline
	void
	projected_parallel( xyzTriple const & v, xyzTriple & a )
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
	distance( xyzTriple const & v ) const
	{
		return std::sqrt( square( x_ - v.x_ ) + square( y_ - v.y_ ) + square( z_ - v.z_ ) );
	}


	/// @brief Distance
	friend
	inline
	Value
	distance( xyzTriple const & a, xyzTriple const & b )
	{
		return std::sqrt( square( a.x_ - b.x_ ) + square( a.y_ - b.y_ ) + square( a.z_ - b.z_ ) );
	}


	/// @brief Distance squared
	inline
	Value
	distance_squared( xyzTriple const & v ) const
	{
		return square( x_ - v.x_ ) + square( y_ - v.y_ ) + square( z_ - v.z_ );
	}


	/// @brief Distance squared
	friend
	inline
	Value
	distance_squared( xyzTriple const & a, xyzTriple const & b )
	{
		return square( a.x_ - b.x_ ) + square( a.y_ - b.y_ ) + square( a.z_ - b.z_ );
	}


	/// @brief Dot product
	inline
	Value
	dot( xyzTriple const & v ) const
	{
		return ( x_ * v.x_ ) + ( y_ * v.y_ ) + ( z_ * v.z_ );
	}


	/// @brief Dot product
	inline
	Value
	dot_product( xyzTriple const & v ) const
	{
		return ( x_ * v.x_ ) + ( y_ * v.y_ ) + ( z_ * v.z_ );
	}


	/// @brief Inner product ( == dot product )
	inline
	Value
	inner_product( xyzTriple const & v ) const
	{
		return ( x_ * v.x_ ) + ( y_ * v.y_ ) + ( z_ * v.z_ );
	}


	/// @brief Dot product
	friend
	inline
	Value
	dot( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.x_ * b.x_ ) + ( a.y_ * b.y_ ) + ( a.z_ * b.z_ );
	}


	/// @brief Dot product
	friend
	inline
	Value
	dot_product( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.x_ * b.x_ ) + ( a.y_ * b.y_ ) + ( a.z_ * b.z_ );
	}


	/// @brief Inner product ( == dot product )
	friend
	inline
	Value
	inner_product( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.x_ * b.x_ ) + ( a.y_ * b.y_ ) + ( a.z_ * b.z_ );
	}


	/// @brief Cross product
	inline
	xyzTriple
	cross( xyzTriple const & v ) const
	{
		return xyzTriple(
			( y_ * v.z_ ) - ( z_ * v.y_ ),
			( z_ * v.x_ ) - ( x_ * v.z_ ),
			( x_ * v.y_ ) - ( y_ * v.x_ )
		);
	}


	/// @brief Cross product
	inline
	xyzTriple
	cross_product( xyzTriple const & v ) const
	{
		return xyzTriple(
			( y_ * v.z_ ) - ( z_ * v.y_ ),
			( z_ * v.x_ ) - ( x_ * v.z_ ),
			( x_ * v.y_ ) - ( y_ * v.x_ )
		);
	}


	/// @brief Cross product
	friend
	inline
	xyzTriple
	cross( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple(
			( a.y_ * b.z_ ) - ( a.z_ * b.y_ ),
			( a.z_ * b.x_ ) - ( a.x_ * b.z_ ),
			( a.x_ * b.y_ ) - ( a.y_ * b.x_ )
		);
	}


	/// @brief Cross product
	friend
	inline
	xyzTriple
	cross_product( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple(
			( a.y_ * b.z_ ) - ( a.z_ * b.y_ ),
			( a.z_ * b.x_ ) - ( a.x_ * b.z_ ),
			( a.x_ * b.y_ ) - ( a.y_ * b.x_ )
		);
	}


	/// @brief Cross product: Return via argument (slightly faster)
	friend
	inline
	void
	cross( xyzTriple const & a, xyzTriple const & b, xyzTriple & c )
	{
		c.x_ = ( a.y_ * b.z_ ) - ( a.z_ * b.y_ );
		c.y_ = ( a.z_ * b.x_ ) - ( a.x_ * b.z_ );
		c.z_ = ( a.x_ * b.y_ ) - ( a.y_ * b.x_ );
	}


	/// @brief Cross product: Return via argument (slightly faster)
	friend
	inline
	void
	cross_product( xyzTriple const & a, xyzTriple const & b, xyzTriple & c )
	{
		c.x_ = ( a.y_ * b.z_ ) - ( a.z_ * b.y_ );
		c.y_ = ( a.z_ * b.x_ ) - ( a.x_ * b.z_ );
		c.z_ = ( a.x_ * b.y_ ) - ( a.y_ * b.x_ );
	}


	/// @brief Midpoint of 2 xyzTriples
	friend
	inline
	xyzTriple
	midpoint( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple(
			Value( 0.5 * ( a.x_ + b.x_ ) ),
			Value( 0.5 * ( a.y_ + b.y_ ) ),
			Value( 0.5 * ( a.z_ + b.z_ ) )
		);
	}


	/// @brief Midpoint of 2 xyzTriples: Return via argument (slightly faster)
	friend
	inline
	void
	midpoint( xyzTriple const & a, xyzTriple const & b, xyzTriple & m )
	{
		m.x_ = Value( 0.5 * ( a.x_ + b.x_ ) );
		m.y_ = Value( 0.5 * ( a.y_ + b.y_ ) );
		m.z_ = Value( 0.5 * ( a.z_ + b.z_ ) );
	}


	/// @brief Center of 2 xyzTriples
	friend
	inline
	xyzTriple
	center( xyzTriple const & a, xyzTriple const & b )
	{
		return xyzTriple(
			Value( 0.5 * ( a.x_ + b.x_ ) ),
			Value( 0.5 * ( a.y_ + b.y_ ) ),
			Value( 0.5 * ( a.z_ + b.z_ ) )
		);
	}


	/// @brief Center of 2 xyzTriples: Return via argument (slightly faster)
	friend
	inline
	void
	center( xyzTriple const & a, xyzTriple const & b, xyzTriple & m )
	{
		m.x_ = Value( 0.5 * ( a.x_ + b.x_ ) );
		m.y_ = Value( 0.5 * ( a.y_ + b.y_ ) );
		m.z_ = Value( 0.5 * ( a.z_ + b.z_ ) );
	}


	/// @brief Center of 3 xyzTriples
	friend
	inline
	xyzTriple
	center( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c )
	{
		long double const third( 1.0 / 3.0 );
		return xyzTriple(
			Value( third * ( a.x_ + b.x_ + c.x_ ) ),
			Value( third * ( a.y_ + b.y_ + c.y_ ) ),
			Value( third * ( a.z_ + b.z_ + c.z_ ) )
		);
	}


	/// @brief Center of 3 xyzTriples: Return via argument (slightly faster)
	friend
	inline
	void
	center( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c, xyzTriple & m )
	{
		long double const third( 1.0 / 3.0 );
		m.x_ = Value( third * ( a.x_ + b.x_ + c.x_ ) );
		m.y_ = Value( third * ( a.y_ + b.y_ + c.y_ ) );
		m.z_ = Value( third * ( a.z_ + b.z_ + c.z_ ) );
	}


	/// @brief Center of 4 xyzTriples
	friend
	inline
	xyzTriple
	center( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c, xyzTriple const & d )
	{
		return xyzTriple(
			Value( 0.25 * ( a.x_ + b.x_ + c.x_ + d.x_ ) ),
			Value( 0.25 * ( a.y_ + b.y_ + c.y_ + d.y_ ) ),
			Value( 0.25 * ( a.z_ + b.z_ + c.z_ + d.z_ ) )
		);
	}


	/// @brief Center of 4 xyzTriples: Return via argument (slightly faster)
	friend
	inline
	void
	center( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c, xyzTriple const & d, xyzTriple & m )
	{
		m.x_ = Value( 0.25 * ( a.x_ + b.x_ + c.x_ + d.x_ ) );
		m.y_ = Value( 0.25 * ( a.y_ + b.y_ + c.y_ + d.y_ ) );
		m.z_ = Value( 0.25 * ( a.z_ + b.z_ + c.z_ + d.z_ ) );
	}


	/// @brief Angle between two vectors (in radians on [ 0, pi ])
	friend
	inline
	Value
	angle_of( xyzTriple const & a, xyzTriple const & b )
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
	angle_of( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c )
	{
		return angle_of( a - b, c - b );
	}


	/// @brief Cosine of angle between two vectors
	friend
	inline
	Value
	cos_of( xyzTriple const & a, xyzTriple const & b )
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
	cos_of( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c )
	{
		return cos_of( a - b, c - b );
	}


	/// @brief Sine of angle between two vectors
	friend
	inline
	Value
	sin_of( xyzTriple const & a, xyzTriple const & b )
	{
		return std::sqrt( Value( 1 ) - square( cos_of( a, b ) ) );
	}


	/// @brief Sine of angle formed by three consecutive points
	/// @note  For points a, b, c, the angle is the angle between the vectors a - b  and c - b
	///        in other words, the positive angle about b from a to c
	friend
	inline
	Value
	sin_of( xyzTriple const & a, xyzTriple const & b, xyzTriple const & c )
	{
		return sin_of( a - b, c - b );
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
			throw std::out_of_range("numeric::xyzTriple::at");
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
			throw std::out_of_range ("numeric::xyzTriple::at");
		}

		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}

	/// @brief xyzTriple[ i ] const: 0-based index
	inline
	Value const &
	operator []( int const i ) const
	{
		assert( ( i >= 0 ) && ( i < 3 ) );
		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}


	/// @brief xyzTriple[ i ]: 0-based index
	inline
	Value &
	operator []( int const i )
	{
		assert( ( i >= 0 ) && ( i < 3 ) );
		return ( i == 0 ? x_ : ( i == 1 ? y_ : z_ ) );
	}


	/// @brief xyzTriple( i ) const: 1-based index
	inline
	Value const &
	operator ()( int const i ) const
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		return ( i == 1 ? x_ : ( i == 2 ? y_ : z_ ) );
	}


	/// @brief xyzTriple( i ): 1-based index
	inline
	Value &
	operator ()( int const i )
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		return ( i == 1 ? x_ : ( i == 2 ? y_ : z_ ) );
	}


public: // Comparison


	/// @brief xyzTriple == xyzTriple
	friend
	inline
	bool
	operator ==( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.x_ == b.x_ ) && ( a.y_ == b.y_ ) && ( a.z_ == b.z_ );
	}


	/// @brief xyzTriple != xyzTriple
	friend
	inline
	bool
	operator !=( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.x_ != b.x_ ) || ( a.y_ != b.y_ ) || ( a.z_ != b.z_ );
	}


	/// @brief xyzTriple < xyzTriple: Lexicographic order
	friend
	inline
	bool
	operator <( xyzTriple const & a, xyzTriple const & b )
	{
		return (
		 ( a.x_ < b.x_ ? true :
		 ( b.x_ < a.x_ ? false : // a.x_ == b.x_
		 ( a.y_ < b.y_ ? true :
		 ( b.y_ < a.y_ ? false : // a.y_ == b.y_
		 ( a.z_ < b.z_ ) ) ) ) ) );
	}


	/// @brief xyzTriple <= xyzTriple
	friend
	inline
	bool
	operator <=( xyzTriple const & a, xyzTriple const & b )
	{
		return (
		 ( a.x_ < b.x_ ? true :
		 ( b.x_ < a.x_ ? false : // a.x_ == b.x_
		 ( a.y_ < b.y_ ? true :
		 ( b.y_ < a.y_ ? false : // a.y_ == b.y_
		 ( a.z_ <= b.z_ ) ) ) ) ) );
	}


	/// @brief xyzTriple >= xyzTriple
	friend
	inline
	bool
	operator >=( xyzTriple const & a, xyzTriple const & b )
	{
		return (
		 ( a.x_ > b.x_ ? true :
		 ( b.x_ > a.x_ ? false : // a.x_ == b.x_
		 ( a.y_ > b.y_ ? true :
		 ( b.y_ > a.y_ ? false : // a.y_ == b.y_
		 ( a.z_ >= b.z_ ) ) ) ) ) );
	}


	/// @brief xyzTriple > xyzTriple
	friend
	inline
	bool
	operator >( xyzTriple const & a, xyzTriple const & b )
	{
		return (
		 ( a.x_ > b.x_ ? true :
		 ( b.x_ > a.x_ ? false : // a.x_ == b.x_
		 ( a.y_ > b.y_ ? true :
		 ( b.y_ > a.y_ ? false : // a.y_ == b.y_
		 ( a.z_ > b.z_ ) ) ) ) ) );
	}


	/// @brief xyzTriple == Value
	friend
	inline
	bool
	operator ==( xyzTriple const & v, Value const & t )
	{
		return ( v.x_ == t ) && ( v.y_ == t ) && ( v.z_ == t );
	}


	/// @brief xyzTriple != Value
	friend
	inline
	bool
	operator !=( xyzTriple const & v, Value const & t )
	{
		return ( v.x_ != t ) || ( v.y_ != t ) || ( v.z_ != t );
	}


	/// @brief xyzTriple < Value
	friend
	inline
	bool
	operator <( xyzTriple const & v, Value const & t )
	{
		return ( v.x_ < t ) && ( v.y_ < t ) && ( v.z_ < t );
	}


	/// @brief xyzTriple <= Value
	friend
	inline
	bool
	operator <=( xyzTriple const & v, Value const & t )
	{
		return ( v.x_ <= t ) && ( v.y_ <= t ) && ( v.z_ <= t );
	}


	/// @brief xyzTriple >= Value
	friend
	inline
	bool
	operator >=( xyzTriple const & v, Value const & t )
	{
		return ( v.x_ >= t ) && ( v.y_ >= t ) && ( v.z_ >= t );
	}


	/// @brief xyzTriple > Value
	friend
	inline
	bool
	operator >( xyzTriple const & v, Value const & t )
	{
		return ( v.x_ > t ) && ( v.y_ > t ) && ( v.z_ > t );
	}


	/// @brief Value == xyzTriple
	friend
	inline
	bool
	operator ==( Value const & t, xyzTriple const & v )
	{
		return ( t == v.x_ ) && ( t == v.y_ ) && ( t == v.z_ );
	}


	/// @brief Value != xyzTriple
	friend
	inline
	bool
	operator !=( Value const & t, xyzTriple const & v )
	{
		return ( t != v.x_ ) || ( t != v.y_ ) || ( t != v.z_ );
	}


	/// @brief Value < xyzTriple
	friend
	inline
	bool
	operator <( Value const & t, xyzTriple const & v )
	{
		return ( t < v.x_ ) && ( t < v.y_ ) && ( t < v.z_ );
	}


	/// @brief Value <= xyzTriple
	friend
	inline
	bool
	operator <=( Value const & t, xyzTriple const & v )
	{
		return ( t <= v.x_ ) && ( t <= v.y_ ) && ( t <= v.z_ );
	}


	/// @brief Value >= xyzTriple
	friend
	inline
	bool
	operator >=( Value const & t, xyzTriple const & v )
	{
		return ( t >= v.x_ ) && ( t >= v.y_ ) && ( t >= v.z_ );
	}


	/// @brief Value > xyzTriple
	friend
	inline
	bool
	operator >( Value const & t, xyzTriple const & v )
	{
		return ( t > v.x_ ) && ( t > v.y_ ) && ( t > v.z_ );
	}


	/// @brief Equal length?
	inline
	bool
	equal_length( xyzTriple const & v )
	{
		return ( length_squared() == v.length_squared() );
	}


	/// @brief Equal length?
	friend
	inline
	bool
	equal_length( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.length_squared() == b.length_squared() );
	}


	/// @brief Not equal length?
	inline
	bool
	not_equal_length( xyzTriple const & v )
	{
		return ( length_squared() != v.length_squared() );
	}


	/// @brief Not equal length?
	friend
	inline
	bool
	not_equal_length( xyzTriple const & a, xyzTriple const & b )
	{
		return ( a.length_squared() != b.length_squared() );
	}


	/// @brief Longer?
	inline
	bool
	longer( xyzTriple const & v )
	{
		return ( length_squared() > v.length_squared() );
	}


	/// @brief Longer or equal length?
	inline
	bool
	longer_or_equal( xyzTriple const & v )
	{
		return ( length_squared() >= v.length_squared() );
	}


	/// @brief Shorter?
	inline
	bool
	shorter( xyzTriple const & v )
	{
		return ( length_squared() < v.length_squared() );
	}


	/// @brief Shorter or equal length?
	inline
	bool
	shorter_or_equal( xyzTriple const & v )
	{
		return ( length_squared() <= v.length_squared() );
	}


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


}; // xyzTriple


/// @brief xyzTriple + xyzTriple
template< typename T >
xyzTriple< T >
operator +( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple + T
template< typename T >
xyzTriple< T >
operator +( xyzTriple< T > const & v, T const & t );


/// @brief T + xyzTriple
template< typename T >
xyzTriple< T >
operator +( T const & t, xyzTriple< T > const & v );


/// @brief xyzTriple - xyzTriple
template< typename T >
xyzTriple< T >
operator -( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple - T
template< typename T >
xyzTriple< T >
operator -( xyzTriple< T > const & v, T const & t );


/// @brief T - xyzTriple
template< typename T >
xyzTriple< T >
operator -( T const & t, xyzTriple< T > const & v );


/// @brief xyzTriple * T
template< typename T >
xyzTriple< T >
operator *( xyzTriple< T > const & v, T const & t );


/// @brief T * xyzTriple
template< typename T >
xyzTriple< T >
operator *( T const & t, xyzTriple< T > const & v );


/// @brief xyzTriple / T
template< typename T >
xyzTriple< T >
operator /( xyzTriple< T > const & v, T const & t );


/// @brief Add: xyzTriple + xyzTriple
template< typename T >
void
add( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > & r );


/// @brief Add: xyzTriple + T
template< typename T >
void
add( xyzTriple< T > const & v, T const & t, xyzTriple< T > & r );


/// @brief Add: T + xyzTriple
template< typename T >
void
add( T const & t, xyzTriple< T > const & v, xyzTriple< T > & r );


/// @brief Subtract: xyzTriple - xyzTriple
template< typename T >
void
subtract( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > & r );


/// @brief Subtract: xyzTriple - T
template< typename T >
void
subtract( xyzTriple< T > const & v, T const & t, xyzTriple< T > & r );


/// @brief Subtract: T - xyzTriple
template< typename T >
void
subtract( T const & t, xyzTriple< T > const & v, xyzTriple< T > & r );


/// @brief Multiply: xyzTriple * T
template< typename T >
void
multiply( xyzTriple< T > const & v, T const & t, xyzTriple< T > & r );


/// @brief Multiply: T * xyzTriple
template< typename T >
void
multiply( T const & t, xyzTriple< T > const & v, xyzTriple< T > & r );


/// @brief Divide: xyzTriple / T
template< typename T >
void
divide( xyzTriple< T > const & v, T const & t, xyzTriple< T > & r );


/// @brief xyzTriple with min coordinates of two xyzTriples
template< typename T >
xyzTriple< T >
min( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple with max coordinates of two xyzTriples
template< typename T >
xyzTriple< T >
max( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Distance
template< typename T >
T
distance( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Distance squared
template< typename T >
T
distance_squared( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Dot product
template< typename T >
T
dot( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Dot product
template< typename T >
T
dot_product( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Inner product ( == dot product )
template< typename T >
T
inner_product( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Cross product
template< typename T >
xyzTriple< T >
cross( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Cross product
template< typename T >
xyzTriple< T >
cross_product( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Cross product: Return via argument (slightly faster)
template< typename T >
void
cross( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > & c );


/// @brief Cross product: Return via argument (slightly faster)
template< typename T >
void
cross_product( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > & c );


/// @brief Midpoint of 2 xyzTriples
template< typename T >
xyzTriple< T >
midpoint( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Midpoint of 2 xyzTriples: Return via argument (slightly faster)
template< typename T >
void
midpoint( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > & m );


/// @brief Center of 2 xyzTriples
template< typename T >
xyzTriple< T >
center( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Center of 2 xyzTriples: Return via argument (slightly faster)
template< typename T >
void
center( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > & m );


/// @brief Center of 3 xyzTriples
template< typename T >
xyzTriple< T >
center( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c );


/// @brief Center of 3 xyzTriples: Return via argument (slightly faster)
template< typename T >
void
center( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c, xyzTriple< T > & m );


/// @brief Center of 4 xyzTriples
template< typename T >
xyzTriple< T >
center( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c, xyzTriple< T > const & d );


/// @brief Center of 4 xyzTriples: Return via argument (slightly faster)
template< typename T >
void
center( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c, xyzTriple< T > const & d, xyzTriple< T > & m );


/// @brief Angle between two vectors (in radians on [ 0, pi ])
template< typename T >
T
angle_of( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Angle formed by three consecutive points (in radians on [ 0, pi ])
template< typename T >
T
angle_of( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c );


/// @brief Cosine of angle between two vectors
template< typename T >
T
cos_of( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Cosine of angle formed by three consecutive points
template< typename T >
T
cos_of( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c );


/// @brief Sine of angle between two vectors
template< typename T >
T
sin_of( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Sine of angle formed by three consecutive points
template< typename T >
T
sin_of( xyzTriple< T > const & a, xyzTriple< T > const & b, xyzTriple< T > const & c );


/// @brief xyzTriple == xyzTriple
template< typename T >
bool
operator ==( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple != xyzTriple
template< typename T >
bool
operator !=( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple < xyzTriple
template< typename T >
bool
operator <( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple <= xyzTriple
template< typename T >
bool
operator <=( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple >= xyzTriple
template< typename T >
bool
operator >=( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple > xyzTriple
template< typename T >
bool
operator >( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief xyzTriple == T
template< typename T >
bool
operator ==( xyzTriple< T > const & v, T const & t );


/// @brief xyzTriple != T
template< typename T >
bool
operator !=( xyzTriple< T > const & v, T const & t );


/// @brief xyzTriple < T
template< typename T >
bool
operator <( xyzTriple< T > const & v, T const & t );


/// @brief xyzTriple <= T
template< typename T >
bool
operator <=( xyzTriple< T > const & v, T const & t );


/// @brief xyzTriple >= T
template< typename T >
bool
operator >=( xyzTriple< T > const & v, T const & t );


/// @brief xyzTriple > T
template< typename T >
bool
operator >( xyzTriple< T > const & v, T const & t );


/// @brief T == xyzTriple
template< typename T >
bool
operator ==( T const & t, xyzTriple< T > const & v );


/// @brief T != xyzTriple
template< typename T >
bool
operator !=( T const & t, xyzTriple< T > const & v );


/// @brief T < xyzTriple
template< typename T >
bool
operator <( T const & t, xyzTriple< T > const & v );


/// @brief T <= xyzTriple
template< typename T >
bool
operator <=( T const & t, xyzTriple< T > const & v );


/// @brief T >= xyzTriple
template< typename T >
bool
operator >=( T const & t, xyzTriple< T > const & v );


/// @brief T > xyzTriple
template< typename T >
bool
operator >( T const & t, xyzTriple< T > const & v );


/// @brief Equal length?
template< typename T >
bool
equal_length( xyzTriple< T > const & a, xyzTriple< T > const & b );


/// @brief Not equal length?
template< typename T >
bool
not_equal_length( xyzTriple< T > const & a, xyzTriple< T > const & b );


// PyRosetta work around for templates classes
//class xyzTriple_Double : public xyzTriple< double >
//{};


} // namespace numeric


#endif // INCLUDED_numeric_xyzTriple_HH
