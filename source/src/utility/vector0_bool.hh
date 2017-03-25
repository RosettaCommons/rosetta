// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector0_bool.hh
/// @brief  vector0: std::vector with assert-checked bounds: bool specialization
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_vector0_bool_HH
#define INCLUDED_utility_vector0_bool_HH


// Unit headers
// Guarantee vector0's presence but avoid circular dependency
#ifndef INCLUDED_utility_vector0_HH
#include <utility/vector0.hh>
#endif // INCLUDED_utility_vector0_HH


namespace utility {


/// @brief std::vector with assert-checked bounds: bool specialization
/// @note
///  @li std::vector with assert-checked bounds in operator[] and a few extras
///  @li Can construct and assign from std::vector and swap with std::vector
///  @li Can compare with std::vector
///  @li Can explicitly convert to std::vector
///  @li Public inheritance from concrete vectorL template is safe here
template< typename A >
class vector0< bool, A > :
	public vectorL< 0, bool, A >
{

private: // Types


	typedef  vectorL< 0, bool, A >  super;


public: // Types


	// STL/boost style
	typedef  typename super::value_type  value_type;
	typedef  typename super::reference  reference;
	typedef  typename super::const_reference  const_reference;
	typedef  typename super::pointer  pointer;
	typedef  typename super::const_pointer  const_pointer;
	typedef  typename super::iterator  iterator;
	typedef  typename super::const_iterator  const_iterator;
	typedef  typename super::reverse_iterator  reverse_iterator;
	typedef  typename super::const_reverse_iterator  const_reverse_iterator;
	typedef  typename super::size_type  size_type;
	typedef  typename super::difference_type  difference_type;
	typedef  typename super::allocator_type  allocator_type;
	typedef  typename super::index_type  index_type;
	typedef  typename super::ssize_type  ssize_type;

	// Project style
	typedef  typename super::Value  Value;
	typedef  typename super::Reference  Reference;
	typedef  typename super::ConstReference  ConstReference;
	typedef  typename super::Pointer  Pointer;
	typedef  typename super::ConstPointer  ConstPointer;
	typedef  typename super::Iterator  Iterator;
	typedef  typename super::ConstIterator  ConstIterator;
	typedef  typename super::ReverseIterator  ReverseIterator;
	typedef  typename super::ConstReverseIterator  ConstReverseIterator;
	typedef  typename super::Size  Size;
	typedef  typename super::Difference  Difference;
	typedef  typename super::Allocator  Allocator;
	typedef  typename super::Index  Index;
	typedef  typename super::SSize  SSize;


public: // Methods: imports


	using super::assign;
	using super::at;
	using super::back;
	using super::begin;
	using super::capacity;
	using super::clear;
	using super::empty;
	using super::end;
	using super::erase;
	using super::flip;
	using super::front;
	using super::get_allocator;
	using super::insert;
	using super::max_size;
	using super::operator [];
	using super::pop_back;
	using super::push_back;
	using super::rbegin;
	using super::rend;
	using super::reserve;
	using super::resize;
	using super::size;
	using super::swap;


public: // Creation


	/// @brief Default constructor
	inline
	explicit
	vector0( allocator_type const & alloc = allocator_type() ) :
		super( alloc )
	{}


	/// @brief Copy constructor
	inline
	vector0( vector0 const & v ) :
		super( v )
	{}


	/// @brief Assignable copy constructor
	template< ssize_type L_, typename T_, typename A_ >
	inline
	vector0( vectorL< L_, T_, A_ > const & v ) :
		super( v.begin(), v.end() )
	{}


	/// @brief std::vector constructor
	inline
	explicit
	vector0( super const & v ) :
		super( v )
	{}


	/// @brief Assignable std::vector constructor
	template< typename T_, typename A_ >
	inline
	explicit
	vector0( std::vector< T_, A_ > const & v ) :
		super( v.begin(), v.end() )
	{}


	/// @brief Size constructor
	inline
	explicit
	vector0( size_type num ) :
		super( num )
	{}


	/// @brief Uniform value constructor
	inline
	vector0(
		size_type const num,
		value_type const & value,
		allocator_type const & alloc = allocator_type()
	) :
		super( num, value, alloc )
	{}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	vector0(
		InputIterator const beg,
		InputIterator const end,
		allocator_type const & alloc = allocator_type()
	) :
		super( beg, end, alloc )
	{}


	/// @brief C++11 Move constructor
	inline
	vector0( vector0 && v ) :
		super( v )
	{}


	/// @brief C++11 vectorL move constructor
	inline
	explicit
	vector0( super && v ) :
		super( v )
	{}

	// We don't fiddle around much with allocators, so I don't bother to do the allocator move constructor specifications

	/// @brief C++11 initializer list constructor
	inline
	vector0( std::initializer_list<bool> init, const A& alloc = A() ) :
		super( init, alloc )
	{}


	/// @brief Destructor
	inline
	virtual
	~vector0()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	vector0 &
	operator =( vector0 const & v )
	{
		if ( this != &v ) {
			super::operator =( v );
		}
		return *this;
	}


	/// @brief Assignable copy assignment
	template< ssize_type L_, typename T_, typename A_ >
	inline
	vector0 &
	operator =( vectorL< L_, T_, A_ > const & v )
	{
		super::assign( v.begin(), v.end() );
		return *this;
	}


	/// @brief std::vector assignment
	inline
	vector0 &
	operator =( super const & v )
	{
		super::operator =( v );
		return *this;
	}


	/// @brief Assignable std::vector assignment
	template< typename T_, typename A_ >
	inline
	vector0 &
	operator =( std::vector< T_, A_ > const & v )
	{
		super::assign( v.begin(), v.end() );
		return *this;
	}


}; // vector0


} // namespace utility


#endif // INCLUDED_utility_vector0_bool_HH
