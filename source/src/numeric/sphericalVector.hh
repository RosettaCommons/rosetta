// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/sphericalVector.hh
/// @brief  Fast spherical-coordinate numeric vector
/// @author Sam DeLuca
///
/// @remarks
///  @li Inline, loop-free functions for speed
///  @li Non-virtual destructor for speed: Not set up for use as a base class
///  @li Pointer constructor and assignment not available for sphericalVectors of pointers
///  @li Numeric vector semantics: spatial partial ordering
///  @li eventually this will have all the functions that xyzVector does, but not today.


#ifndef INCLUDED_numeric_sphericalVector_hh
#define INCLUDED_numeric_sphericalVector_hh


// Unit headers
#include <numeric/sphericalVector.fwd.hh>
// Package headers
#include <numeric/constants.hh>
// C++ headers
#include <utility/assert.hh>
#include <cmath>


namespace numeric {


/// @brief sphericalVector: Fast spherical-coordinate numeric vector
template< typename T >
class sphericalVector
{


private: // Friends

	// Friend functions (for speed of non-inlining debug builds)


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
	sphericalVector()
	{}


	/// @brief Copy constructor
	inline
	sphericalVector( sphericalVector const & v ) :
		phi_( v.phi_ ),
		theta_( v.theta_ ),
		radius_( v.radius_ )
	{}


	/// @brief Copy constructor
	template< typename U >
	inline
	sphericalVector( sphericalVector< U > const & v ) :
		phi_( v.phi_ ),
		theta_( v.theta_ ),
		radius_( v.radius_ )
	{}

	/// @brief Triple value constructor
	inline
	sphericalVector(
		Value const & phi_a,
		Value const & theta_a,
		Value const & radius_a
	 ) :
		phi_( phi_a ),
		theta_( theta_a ),
		radius_( radius_a )
	{}


	/// @brief Pointer to contiguous values constructor
	/// @note  U must be assignable to a Value
	/// @warning No way to check that argument points to three values
	/// @warning Argument missing an & operator will quietly call the uniform value constructor
	template< typename U >
	inline
	explicit
	sphericalVector( U const * p ) :
		phi_( *p ),
		theta_( *++p ),
		radius_( *++p )
	{}


	/// @brief Destructor
	inline
	~sphericalVector()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	sphericalVector &
	operator =( sphericalVector const & v )
	{
		if ( this != &v ) {
			phi_ = v.phi_;
			theta_ = v.theta_;
			radius_ = v.radius_;
		}
		return *this;
	}


	/// @brief Copy assignment
	template< typename U >
	inline
	sphericalVector &
	operator =( sphericalVector< U > const & v )
	{
		phi_ = v.phi_;
		theta_ = v.theta_;
		radius_ = v.radius_;
		return *this;
	}


	/// @brief Assignment from pointer to contiguous values
	/// @warning No way to check that argument points to three values
	template< typename U >
	inline
	sphericalVector &
	operator =( U const * p )
	{
		phi_ = *p;
		theta_ = *++p;
		radius_ = *++p;
		return *this;
	}


	/// @brief += xyzVector
	template< typename U >
	inline
	sphericalVector &
	operator +=( sphericalVector< U > const & v )
	{
		phi_ += v.phi_;
		theta_ += v.theta_;
		radius_ += v.radius_;
		return *this;
	}


	/// @brief -= xyzVector
	template< typename U >
	inline
	sphericalVector &
	operator -=( sphericalVector< U > const & v )
	{
		phi_ -= v.phi_;
		theta_ -= v.theta_;
		radius_ -= v.radius_;
		return *this;
	}


	/// @brief Assign Value * xyzVector
	/// @note  Avoids temporary of = Value * xyzVector
	template< typename U >
	inline
	sphericalVector &
	scaled_assign( Value const & t, sphericalVector< U > const & v )
	{
		phi_ = t * v.phi_;
		theta_ = t * v.theta_;
		radius_ = t * v.radius_;
		return *this;
	}


	/// @brief Add Value * xyzVector
	/// @note  Avoids temporary of += Value * xyzVector
	template< typename U >
	inline
	sphericalVector &
	scaled_add( Value const & t, sphericalVector< U > const & v )
	{
		phi_ += t * v.phi_;
		theta_ += t * v.theta_;
		radius_ += t * v.radius_;
		return *this;
	}


	/// @brief Subtract Value * xyzVector
	/// @note  Avoids temporary of -= Value * xyzVector
	template< typename U >
	inline
	sphericalVector &
	scaled_sub( Value const & t, sphericalVector< U > const & v )
	{
		phi_ -= t * v.phi_;
		theta_ -= t * v.theta_;
		radius_ -= t * v.radius_;
		return *this;
	}


	/// @brief = Value
	inline
	sphericalVector &
	operator =( Value const & t )
	{
		phi_ = theta_ = radius_ = t;
		return *this;
	}


	/// @brief += Value
	inline
	sphericalVector &
	operator +=( Value const & t )
	{
		phi_ += t;
		theta_ += t;
		radius_ += t;
		return *this;
	}


	/// @brief -= Value
	inline
	sphericalVector &
	operator -=( Value const & t )
	{
		phi_ -= t;
		theta_ -= t;
		radius_ -= t;
		return *this;
	}


	/// @brief *= Value
	inline
	sphericalVector &
	operator *=( Value const & t )
	{
		phi_ *= t;
		theta_ *= t;
		radius_ *= t;
		return *this;
	}


	/// @brief /= Value
	inline
	sphericalVector &
	operator /=( Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		phi_ *= inv_t;
		theta_ *= inv_t;
		radius_ *= inv_t;
		return *this;
	}


	/// @brief Triple value assignment
	inline
	sphericalVector &
	assign(
		Value const & phi_a,
		Value const & theta_a,
		Value const & radius_a
	)
	{
		phi_ = phi_a;
		theta_ = theta_a;
		radius_ = radius_a;
		return *this;
	}


public: // Methods


	/// @brief Clear
	inline
	sphericalVector &
	clear()
	{
		phi_ = theta_ = radius_ = Value( 0 );
		return *this;
	}


	/// @brief Zero
	inline
	sphericalVector &
	zero()
	{
		phi_ = theta_ = radius_ = Value( 0 );
		return *this;
	}


	/// @brief sphericalVector + sphericalVector
	friend
	inline
	sphericalVector
	operator +( sphericalVector const & a, sphericalVector const & b )
	{
		return sphericalVector( a.phi_ + b.phi_, a.theta_ + b.theta_, a.radius_ + b.radius_ );
	}


	/// @brief sphericalVector + Value
	friend
	inline
	sphericalVector
	operator +( sphericalVector const & v, Value const & t )
	{
		return sphericalVector( v.phi_ + t, v.theta_ + t, v.radius_ + t );
	}


	/// @brief Value + sphericalVector
	friend
	inline
	sphericalVector
	operator +( Value const & t, sphericalVector const & v )
	{
		return sphericalVector( t + v.phi_, t + v.theta_, t + v.radius_ );
	}


	/// @brief sphericalVector - sphericalVector
	friend
	inline
	sphericalVector
	operator -( sphericalVector const & a, sphericalVector const & b )
	{
		return sphericalVector( a.phi_ - b.phi_, a.theta_ - b.theta_, a.radius_ - b.radius_ );
	}


	/// @brief sphericalVector - Value
	friend
	inline
	sphericalVector
	operator -( sphericalVector const & v, Value const & t )
	{
		return sphericalVector( v.phi_ - t, v.theta_ - t, v.radius_ - t );
	}


	/// @brief Value - sphericalVector
	friend
	inline
	sphericalVector
	operator -( Value const & t, sphericalVector const & v )
	{
		return sphericalVector( t - v.phi_, t - v.theta_, t - v.radius_ );
	}


	/// @brief sphericalVector * Value
	friend
	inline
	sphericalVector
	operator *( sphericalVector const & v, Value const & t )
	{
		return sphericalVector( v.phi_ * t, v.theta_ * t, v.radius_ * t );
	}


	/// @brief Value * xyzVector
	friend
	inline
	sphericalVector
	operator *( Value const & t, sphericalVector const & v )
	{
		return sphericalVector( t * v.phi_, t * v.theta_, t * v.radius_ );
	}


	/// @brief xyzVector / Value
	friend
	inline
	sphericalVector
	operator /( sphericalVector const & v, Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value ( 1 ) / t );
		return sphericalVector( v.phi_ * inv_t, v.theta_ * inv_t, v.radius_ * inv_t );
	}


	/// @brief Add: xyzVector + xyzVector
	friend
	inline
	void
	add( sphericalVector const & a, sphericalVector const & b, sphericalVector & r )
	{
		r.phi_ = a.phi_ + b.phi_;
		r.theta_ = a.theta_ + b.theta_;
		r.radius_ = a.radius_ + b.radius_;
	}


	/// @brief Add: xyzVector + Value
	friend
	inline
	void
	add( sphericalVector const & v, Value const & t, sphericalVector & r )
	{
		r.phi_ = v.phi_ + t;
		r.theta_ = v.theta_ + t;
		r.radius_ = v.radius_ + t;
	}


	/// @brief Add: Value + xyzVector
	friend
	inline
	void
	add( Value const & t, sphericalVector const & v, sphericalVector & r )
	{
		r.phi_ = t + v.phi_;
		r.theta_ = t + v.theta_;
		r.radius_ = t + v.radius_;
	}


	/// @brief Subtract: xyzVector - xyzVector
	friend
	inline
	void
	subtract( sphericalVector const & a, sphericalVector const & b, sphericalVector & r )
	{
		r.phi_ = a.phi_ - b.phi_;
		r.theta_ = a.theta_ - b.theta_;
		r.radius_ = a.radius_ - b.radius_;
	}


	/// @brief Subtract: xyzVector - Value
	friend
	inline
	void
	subtract( sphericalVector const & v, Value const & t, sphericalVector & r )
	{
		r.phi_ = v.phi_ - t;
		r.theta_ = v.theta_ - t;
		r.radius_ = v.radius_ - t;
	}


	/// @brief Subtract: Value - sphericalVector
	friend
	inline
	void
	subtract( Value const & t, sphericalVector const & v, sphericalVector & r )
	{
		r.phi_ = t - v.phi_;
		r.theta_ = t - v.theta_;
		r.radius_ = t - v.radius_;
	}


	/// @brief Multiply: xyzVector * Value
	friend
	inline
	void
	multiply( sphericalVector const & v, Value const & t, sphericalVector & r )
	{
		r.phi_ = v.phi_ * t;
		r.theta_ = v.theta_ * t;
		r.radius_ = v.radius_ * t;
	}


	/// @brief Multiply: Value * xyzVector
	friend
	inline
	void
	multiply( Value const & t, sphericalVector const & v, sphericalVector & r )
	{
		r.phi_ = t * v.phi_;
		r.theta_ = t * v.theta_;
		r.radius_ = t * v.radius_;
	}


	/// @brief Divide: xyzVector / Value
	friend
	inline
	void
	divide( sphericalVector const & v, Value const & t, sphericalVector & r )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		r.phi_ = v.phi_ * inv_t;
		r.theta_ = v.theta_ * inv_t;
		r.radius_ = v.radius_ * inv_t;
	}


	/// @brief Set minimum coordinates wrt another xyzVector
	inline
	sphericalVector &
	min( sphericalVector const & v )
	{
		phi_ = ( phi_ <= v.phi_ ? phi_ : v.phi_ );
		theta_ = ( theta_ <= v.theta_ ? theta_ : v.theta_ );
		radius_ = ( radius_ <= v.radius_ ? radius_ : v.radius_ );
		return *this;
	}


	/// @brief Set maximum coordinates wrt another xyzVector
	inline
	sphericalVector &
	max( sphericalVector const & v )
	{
		phi_ = ( phi_ >= v.phi_ ? phi_ : v.phi_ );
		theta_ = ( theta_ >= v.theta_ ? theta_ : v.theta_ );
		radius_ = ( radius_ >= v.radius_ ? radius_ : v.radius_ );
		return *this;
	}


	/// @brief xyzVector with min coordinates of two xyzVectors
	friend
	inline
	sphericalVector
	min( sphericalVector const & a, sphericalVector const & b )
	{
		return sphericalVector(
		 ( a.phi_ <= b.phi_ ? a.phi_ : b.phi_ ),
		 ( a.theta_ <= b.theta_ ? a.theta_ : b.theta_ ),
		 ( a.radius_ <= b.radius_ ? a.radius_ : b.radius_ )
		);
	}


	/// @brief sphericalVector with max coordinates of two sphericalVector
	friend
	inline
	sphericalVector
	max( sphericalVector const & a, sphericalVector const & b )
	{
		return sphericalVector(
		 ( a.phi_ >= b.phi_ ? a.phi_ : b.phi_ ),
		 ( a.theta_ >= b.theta_ ? a.theta_ : b.theta_ ),
		 ( a.radius_ >= b.radius_ ? a.radius_ : b.radius_ )
		);
	}


public: // Properties: accessors


	/// @brief Value x const
	inline
	Value const &
	phi() const
	{
		return phi_;
	}


	/// @brief Value x
	inline
	Value &
	phi()
	{
		return phi_;
	}


	/// @brief Value y const
	inline
	Value const &
	theta() const
	{
		return theta_;
	}


	/// @brief Value y
	inline
	Value &
	theta()
	{
		return theta_;
	}


	/// @brief Value z const
	inline
	Value const &
	radius() const
	{
		return radius_;
	}


	/// @brief Value z
	inline
	Value &
	radius()
	{
		return radius_;
	}


public: // Properties: value assignment


	/// @brief x assignment
	inline
	void
	phi( Value const & phi_a )
	{
		phi_ = phi_a;
	}


	/// @brief y assignment
	inline
	void
	theta( Value const & theta_a )
	{
		theta_ = theta_a;
	}


	/// @brief z assignment
	inline
	void
	radius( Value const & radius_a )
	{
		radius_ = radius_a;
	}


public: // Comparison


	/// @brief sphericalVector == sphericalVector
	friend
	inline
	bool
	operator ==( sphericalVector const & a, sphericalVector const & b )
	{
		return ( a.phi_ == b.phi_ ) && ( a.theta_ == b.theta_ ) && ( a.radius_ == b.radius_ );
	}


	/// @brief sphericalVector != sphericalVector
	friend
	inline
	bool
	operator !=( sphericalVector const & a, sphericalVector const & b )
	{
		return ( a.phi_ != b.phi_ ) || ( a.theta_ != b.theta_ ) || ( a.radius_ != b.radius_ );
	}


private: // Fields


	/// @brief Coordinates of the 3 coordinate vector
	Value phi_;
	Value theta_;
	Value radius_;


}; // sphericalVector


} // namespace numeric


#endif // INCLUDED_numeric_sphericalVector_HH
