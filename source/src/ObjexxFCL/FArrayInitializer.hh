#ifndef INCLUDED_ObjexxFCL_FArrayInitializer_hh
#define INCLUDED_ObjexxFCL_FArrayInitializer_hh


// FArrayInitializer: FArray Initializer Class Template
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>

// C++ Headers
#include <cassert>


namespace ObjexxFCL {


/// @brief FArrayInitializer: FArray Initializer Class Template
template<
	typename T,
	template< typename > class Array
>
class FArrayInitializer
{


public: // Types


	typedef  FArrayTraits< T >  Traits;

	// STL style
	typedef  T  value_type;
	typedef  void (*function_type)( Array< T > & );

	// C++ style
	typedef  T  Value;
	typedef  void (*Function)( Array< T > & );


private: // Types


	enum State { INACTIVE, VALUE, FUNCTION };


public: // Creation


	/// @brief Default Constructor
	inline
	FArrayInitializer() :
		state_( INACTIVE ),
		value_( Traits::initial_value() ),
		function_( static_cast< function_type >( 0 ) )
	{}


	/// @brief Value Constructor
	inline
	explicit
	FArrayInitializer( T const & value_a ) :
		state_( VALUE ),
		value_( value_a ),
		function_( static_cast< function_type >( 0 ) )
	{}


	/// @brief Function Constructor
	inline
	explicit
	FArrayInitializer( function_type const & function_a ) :
		state_( function_a ? FUNCTION : INACTIVE ),
		value_( Traits::initial_value() ),
		function_( function_a ? function_a : static_cast< function_type >( 0 ) )
	{}


public: // Assignment


	/// @brief Value Assignment
	inline
	FArrayInitializer &
	operator =( T const & value_a )
	{
		state_ = VALUE;
		value_ = value_a;
		function_ = static_cast< function_type >( 0 );
		return *this;
	}


	/// @brief Function Assignment
	inline
	FArrayInitializer &
	operator =( function_type const & function_a )
	{
		state_ = ( function_a ? FUNCTION : INACTIVE );
		value_ = Traits::initial_value();
		function_ = ( function_a ? function_a : static_cast< function_type >( 0 ) );
		return *this;
	}


public: // Inspector


	/// @brief Active?
	inline
	bool
	is_active() const
	{
		return ( state_ != INACTIVE );
	}


	/// @brief Value?
	inline
	bool
	is_value() const
	{
		return ( state_ == VALUE );
	}


	/// @brief Function?
	inline
	bool
	is_function() const
	{
		return ( state_ == FUNCTION );
	}


	/// @brief Value
	inline
	T const &
	value() const
	{
		assert( state_ == VALUE );
		return value_;
	}


	/// @brief Function
	inline
	function_type const &
	function() const
	{
		assert( state_ == FUNCTION );
		return function_;
	}


public: // Modifier


	/// @brief Clear
	inline
	void
	clear()
	{
		state_ = INACTIVE;
		value_ = Traits::initial_value();
		function_ = static_cast< function_type >( 0 );
	}


private: // Data


	/// @brief State
	State state_;

	/// @brief Value
	T value_;

	/// @brief Function
	function_type function_;


}; // FArrayInitializer


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArrayInitializer_HH
