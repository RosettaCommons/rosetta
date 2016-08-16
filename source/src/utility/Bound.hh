// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/Bound.hh
/// @brief  Bound value class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_Bound_hh
#define INCLUDED_utility_Bound_hh


// Unit headers
#include <utility/Bound.fwd.hh>

// C++ headers
#include <utility/assert.hh>


namespace utility {


/// @brief Bound value class
template< typename T >
class Bound
{


public: // Types


	// Project style
	typedef  T  Value;

	// STL/Boost style
	typedef  T  value_type;


public: // Creation


	/// @brief Default constructor
	inline
	Bound() :
		active_( false ),
		value_( Value() ),
		strict_( false )
	{}


	// Value constructor
	inline
	explicit
	Bound(
		Value const & value_a,
		bool const strict_a = false
	) :
		active_( true ),
		value_( value_a ),
		strict_( strict_a )
	{}


	// Strict named constructor
	inline
	static
	Bound
	Strict( Value const & value_a )
	{
		return Bound( value_a, true );
	}


	/// @brief Destructor
	inline
	~Bound()
	{}


public: // Methods


	/// @brief Value assignment
	inline
	Bound &
	value(
		Value const & value_a,
		bool const strict_a = false
	)
	{
		active_ = true;
		value_ = value_a;
		strict_ = strict_a;
		return *this;
	}


	/// @brief Value assignment
	inline
	Bound &
	operator ()(
		Value const & value_a,
		bool const strict_a = false
	)
	{
		active_ = true;
		value_ = value_a;
		strict_ = strict_a;
		return *this;
	}


	/// @brief Activate
	inline
	Bound &
	activate()
	{
		active_ = true;
		return *this;
	}


	/// @brief Deactivate
	inline
	Bound &
	deactivate()
	{
		active_ = false;
		return *this;
	}


	/// @brief Clear
	inline
	Bound &
	clear()
	{
		active_ = false;
		value_ = Value();
		strict_ = false;
		return *this;
	}


public: // Properties


	/// @brief Active bound?
	inline
	bool
	active() const
	{
		return active_;
	}


	/// @brief Inactive bound?
	inline
	bool
	inactive() const
	{
		return ( ! active_ );
	}


	/// @brief Bound value
	inline
	Value const &
	operator ()() const
	{
		debug_assert( active_ );
		return value_;
	}


	/// @brief Bound value
	inline
	Value const &
	value() const
	{
		debug_assert( active_ );
		return value_;
	}


	/// @brief Strict inequality (< or >) bound?
	inline
	bool
	strict() const
	{
		return strict_;
	}


private: // Fields


	/// @brief Active bound?
	bool active_;

	/// @brief Bound value
	Value value_;

	/// @brief Strict inequality (< or >) bound?
	bool strict_;


}; // Bound


} // namespace utility


#endif // INCLUDED_utility_Bound_HH
