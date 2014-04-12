// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vectorL.hh
/// @brief  vectorL: std::vector with L-based indexing
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_vectorL_hh
#define INCLUDED_utility_vectorL_hh


// Unit headers
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL_Selector.hh>

// C++ headers
#include <cassert>
#include <vector>
#include <algorithm>
#include <memory>


namespace utility {


/// @brief std::vector with L-based indexing
/// @note
///  @li std::vector with L-based indexing and a few extras
///  @li Lower index must be in the range of ssize_t
///  @li Index type is std::size_t or ssize_t depending on sign of L
///  @li When L is negative indexing operators can only reach the first max( ssize_t )
///      element and attempting to index beyond that will trigger an assertion failure
///  @li Can construct and assign from std::vector and swap with std::vector
///  @li Can compare with std::vector: compares contents ignoring indexes
///  @li Can explicitly convert to std::vector
///  @li Private inheritance from std::vector is safe here
template< platform::SSize L, typename T, typename A >
class vectorL :
	private std::vector< T, A >
{


private: // Types


	typedef  std::vector< T, A >  super;


protected: // Types


	typedef  std::vector< T, A >  root;


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
	typedef  typename vectorL_IndexSelector< L >= 0 >::index_type  index_type;
	typedef  platform::SSize  ssize_type;

	// Project style
	typedef  typename super::value_type  Value;
	typedef  typename super::reference  Reference;
	typedef  typename super::const_reference  ConstReference;
	typedef  typename super::pointer  Pointer;
	typedef  typename super::const_pointer  ConstPointer;
	typedef  typename super::iterator  Iterator;
	typedef  typename super::const_iterator  ConstIterator;
	typedef  typename super::reverse_iterator  ReverseIterator;
	typedef  typename super::const_reverse_iterator  ConstReverseIterator;
	typedef  typename super::size_type  Size;
	typedef  typename super::difference_type  Difference;
	typedef  typename super::allocator_type  Allocator;
	typedef  typename vectorL_IndexSelector< L >= 0 >::Index  Index;
	typedef  platform::SSize  SSize;


public: // Methods: imports


	using super::assign;
	using super::back;
	using super::begin;
	using super::capacity;
	using super::clear;
	using super::empty;
	using super::end;
	using super::erase;
	using super::front;
	using super::get_allocator;
	using super::insert;
	using super::max_size;
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
	vectorL( allocator_type const & alloc = allocator_type() ) :
		super( alloc )
	{}


	/// @brief Copy constructor
	inline
	vectorL( vectorL const & v ) :
		super( v )
	{}


	/// @brief Assignable copy constructor
	template< ssize_type L_, typename T_, typename A_ >
	inline
	vectorL( vectorL< L_, T_, A_ > const & v ) :
		super( v.begin(), v.end() )
	{}


	/// @brief std::vector constructor
	inline
	explicit
	vectorL( super const & v ) :
		super( v )
	{}


	/// @brief Assignable std::vector constructor
	template< typename T_, typename A_ >
	inline
	explicit
	vectorL( std::vector< T_, A_ > const & v ) :
		super( v.begin(), v.end() )
	{}


	/// @brief Size constructor
	inline
	explicit
	vectorL( size_type const num ) :
		super( num )
	{}


	/// @brief Uniform value constructor
	inline
	vectorL(
		size_type const num,
		value_type const & value,
		allocator_type const & alloc = allocator_type()
	) :
		super( num, value, alloc )
	{}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	vectorL(
		InputIterator const beg,
		InputIterator const end,
		allocator_type const & alloc = allocator_type()
	) :
		super( beg, end, alloc )
	{}


	/// @brief Destructor
	inline
	virtual
	~vectorL()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	vectorL &
	operator =( vectorL const & v )
	{
		if ( this != &v ) {
			super::operator =( v );
		}
		return *this;
	}


	/// @brief Assignable copy assignment
	template< ssize_type L_, typename T_, typename A_ >
	inline
	vectorL &
	operator =( vectorL< L_, T_, A_ > const & v )
	{
		super::assign( v.begin(), v.end() );
		return *this;
	}


	/// @brief std::vector assignment
	inline
	vectorL &
	operator =( super const & v )
	{
		super::operator =( v );
		return *this;
	}


	/// @brief Assignable std::vector assignment
	template< typename T_, typename A_ >
	inline
	vectorL &
	operator =( std::vector< T_, A_ > const & v )
	{
		super::assign( v.begin(), v.end() );
		return *this;
	}


public: // Conversion


	/// @brief std::vector const explicit conversion
	inline
	super const &
	vector() const
	{
		return static_cast< super const & >( *this );
	}


	/// @brief std::vector explicit conversion
	inline
	super &
	vector()
	{
		return static_cast< super & >( *this );
	}


public: // Methods


	/// @brief Add an element to the back of the vector
	inline
	vectorL &
	add_back( T const & t )
	{
		push_back( t );
		return *this;
	}


	/// @brief Remove the element at the back of the vector
	inline
	vectorL &
	remove_back()
	{
		pop_back();
		return *this;
	}


	/// @brief Shrink the index map to remove unused capacity
	inline
	void
	shrink()
	{
		if ( super::size() < super::capacity() ) vectorL( *this ).swap( *this );
	}


	/// @brief  Check if vector contains a given element.
	inline
	bool
	contains(T const & t) const
	{
		if (std::find(begin(), end(), t) == end()) {
			return false;
		}
		else {
			return true;
		}
	}


	/// @brief  Return the index of a given element or NULL if not found.
	inline
	index_type
	index_of(T const & t)
	{
		index_type loc = std::find(begin(), end(), t) - begin();
		if (loc < size()) {
			return loc + l_;
		}
		else {
			return NULL;
		}
	}

public: // Properties


	/// @brief Has an element with an index?
	inline
	bool
	has( index_type const i ) const
	{
		return ( ( vectorL_ZeroSelector< L != 0 >::ge( i, l_ ) ) && ( static_cast< size_type >( i - l_ ) < super::size() ) );
	}


public: // Indexers


	/// @brief vectorL[ i ] const
	inline
	const_reference
	operator []( index_type const i ) const
	{
		assert( vectorL_ZeroSelector< L != 0 >::ge( i, l_ ) ); // Avoid "always true" warnings when L==0
		assert( static_cast< size_type >( i - l_ ) < super::size() ); // Upper bound check
		return super::operator []( i - l_ );
	}


	/// @brief vectorL[ i ]
	inline
	reference
	operator []( index_type const i )
	{
		assert( vectorL_ZeroSelector< L != 0 >::ge( i, l_ ) ); // Avoid "always true" warnings when L==0
		assert( static_cast< size_type >( i - l_ ) < super::size() ); // Upper bound check
		return super::operator []( i - l_ );
	}


	/// @brief vectorL.at( i ) const
	inline
	const_reference
	at( index_type const i ) const
	{
		assert( vectorL_ZeroSelector< L != 0 >::ge( i, l_ ) ); // Avoid "always true" warnings when L==0
		return super::at( i - l_ );
	}


	/// @brief vectorL.at( i )
	inline
	reference
	at( index_type const i )
	{
		assert( vectorL_ZeroSelector< L != 0 >::ge( i, l_ ) ); // Avoid "always true" warnings when L==0
		return super::at( i - l_ );
	}


	/// @brief Lower index
	inline
	index_type
	l() const
	{
		return l_;
	}


	/// @brief Upper index
	inline
	index_type
	u() const
	{
		assert( ! super::empty() ); // Upper index only meaningful for non-empty vectors
		assert( static_cast< index_type >( super::size() ) >= 0 ); // Catch size range error
		return l_ + static_cast< index_type >( super::size() ) - 1;
	}


public: // Comparison


	/// @brief vectorL == vectorL
	friend
	inline
	bool
	operator ==( vectorL const & a, vectorL const & b )
	{
		return ( static_cast< super const & >( a ) == static_cast< super const & >( b ) );
	}


	/// @brief vectorL != vectorL
	friend
	inline
	bool
	operator !=( vectorL const & a, vectorL const & b )
	{
		return ( static_cast< super const & >( a ) != static_cast< super const & >( b ) );
	}


	/// @brief vectorL < vectorL
	friend
	inline
	bool
	operator <( vectorL const & a, vectorL const & b )
	{
		return ( static_cast< super const & >( a ) < static_cast< super const & >( b ) );
	}


	/// @brief vectorL <= vectorL
	friend
	inline
	bool
	operator <=( vectorL const & a, vectorL const & b )
	{
		return ( static_cast< super const & >( a ) <= static_cast< super const & >( b ) );
	}


	/// @brief vectorL >= vectorL
	friend
	inline
	bool
	operator >=( vectorL const & a, vectorL const & b )
	{
		return ( static_cast< super const & >( a ) >= static_cast< super const & >( b ) );
	}


	/// @brief vectorL > vectorL
	friend
	inline
	bool
	operator >( vectorL const & a, vectorL const & b )
	{
		return ( static_cast< super const & >( a ) > static_cast< super const & >( b ) );
	}


	/// @brief vectorL == std::vector
	friend
	inline
	bool
	operator ==( vectorL const & a, super const & b )
	{
		return ( static_cast< super const & >( a ) == b );
	}


	/// @brief vectorL != std::vector
	friend
	inline
	bool
	operator !=( vectorL const & a, super const & b )
	{
		return ( static_cast< super const & >( a ) != b );
	}


	/// @brief vectorL < std::vector
	friend
	inline
	bool
	operator <( vectorL const & a, super const & b )
	{
		return ( static_cast< super const & >( a ) < b );
	}


	/// @brief vectorL <= std::vector
	friend
	inline
	bool
	operator <=( vectorL const & a, super const & b )
	{
		return ( static_cast< super const & >( a ) <= b );
	}


	/// @brief vectorL >= std::vector
	friend
	inline
	bool
	operator >=( vectorL const & a, super const & b )
	{
		return ( static_cast< super const & >( a ) >= b );
	}


	/// @brief vectorL > std::vector
	friend
	inline
	bool
	operator >( vectorL const & a, super const & b )
	{
		return ( static_cast< super const & >( a ) > b );
	}


	/// @brief std::vector == vectorL
	friend
	inline
	bool
	operator ==( super const & a, vectorL const & b )
	{
		return ( a == static_cast< super const & >( b ) );
	}


	/// @brief std::vector != vectorL
	friend
	inline
	bool
	operator !=( super const & a, vectorL const & b )
	{
		return ( a != static_cast< super const & >( b ) );
	}


	/// @brief std::vector < vectorL
	friend
	inline
	bool
	operator <( super const & a, vectorL const & b )
	{
		return ( a < static_cast< super const & >( b ) );
	}


	/// @brief std::vector <= vectorL
	friend
	inline
	bool
	operator <=( super const & a, vectorL const & b )
	{
		return ( a <= static_cast< super const & >( b ) );
	}


	/// @brief std::vector >= vectorL
	friend
	inline
	bool
	operator >=( super const & a, vectorL const & b )
	{
		return ( a >= static_cast< super const & >( b ) );
	}


	/// @brief std::vector > vectorL
	friend
	inline
	bool
	operator >( super const & a, vectorL const & b )
	{
		return ( a > static_cast< super const & >( b ) );
	}


public: // Swap


	/// @brief swap( vectorL )
	inline
	void
	swap( vectorL & v )
	{
		super::swap( static_cast< super & >( v ) );
	}


	/// @brief swap( vectorL, vectorL )
	friend
	inline
	void
	swap( vectorL & a, vectorL & b )
	{
		static_cast< super & >( a ).swap( static_cast< super & >( b ) );
	}


	/// @brief swap( vectorL, std::vector )
	friend
	inline
	void
	swap( vectorL & a, super & b )
	{
		static_cast< super & >( a ).swap( b );
	}


	/// @brief swap( std::vector, vectorL )
	friend
	inline
	void
	swap( super & a, vectorL & b )
	{
		a.swap( static_cast< super & >( b ) );
	}


private: // Static fields


	/// @brief Lower index in index type
	static index_type const l_;


}; // vectorL


// Static field definitions
template< platform::SSize L, typename T, typename A >
typename vectorL< L, T, A >::index_type const vectorL< L, T, A >::l_( L );


// Friend function namespace declarations


/// @brief vectorL == vectorL
template< platform::SSize L, typename T, typename A >
bool
operator ==( vectorL< L, T, A > const & a, vectorL< L, T, A > const & b );


/// @brief vectorL != vectorL
template< platform::SSize L, typename T, typename A >
bool
operator !=( vectorL< L, T, A > const & a, vectorL< L, T, A > const & b );


/// @brief vectorL < vectorL
template< platform::SSize L, typename T, typename A >
bool
operator <( vectorL< L, T, A > const & a, vectorL< L, T, A > const & b );


/// @brief vectorL <= vectorL
template< platform::SSize L, typename T, typename A >
bool
operator <=( vectorL< L, T, A > const & a, vectorL< L, T, A > const & b );


/// @brief vectorL >= vectorL
template< platform::SSize L, typename T, typename A >
bool
operator >=( vectorL< L, T, A > const & a, vectorL< L, T, A > const & b );


/// @brief vectorL > vectorL
template< platform::SSize L, typename T, typename A >
bool
operator >( vectorL< L, T, A > const & a, vectorL< L, T, A > const & b );


/// @brief vectorL == std::vector
template< platform::SSize L, typename T, typename A >
bool
operator ==( vectorL< L, T, A > const & a, std::vector< T, A > const & b );


/// @brief vectorL != std::vector
template< platform::SSize L, typename T, typename A >
bool
operator !=( vectorL< L, T, A > const & a, std::vector< T, A > const & b );


/// @brief vectorL < std::vector
template< platform::SSize L, typename T, typename A >
bool
operator <( vectorL< L, T, A > const & a, std::vector< T, A > const & b );


/// @brief vectorL <= std::vector
template< platform::SSize L, typename T, typename A >
bool
operator <=( vectorL< L, T, A > const & a, std::vector< T, A > const & b );


/// @brief vectorL >= std::vector
template< platform::SSize L, typename T, typename A >
bool
operator >=( vectorL< L, T, A > const & a, std::vector< T, A > const & b );


/// @brief vectorL > std::vector
template< platform::SSize L, typename T, typename A >
bool
operator >( vectorL< L, T, A > const & a, std::vector< T, A > const & b );


/// @brief std::vector == vectorL
template< platform::SSize L, typename T, typename A >
bool
operator ==( std::vector< T, A > const & a, vectorL< L, T, A > const & b );


/// @brief std::vector != vectorL
template< platform::SSize L, typename T, typename A >
bool
operator !=( std::vector< T, A > const & a, vectorL< L, T, A > const & b );


/// @brief std::vector < vectorL
template< platform::SSize L, typename T, typename A >
bool
operator <( std::vector< T, A > const & a, vectorL< L, T, A > const & b );


/// @brief std::vector <= vectorL
template< platform::SSize L, typename T, typename A >
bool
operator <=( std::vector< T, A > const & a, vectorL< L, T, A > const & b );


/// @brief std::vector >= vectorL
template< platform::SSize L, typename T, typename A >
bool
operator >=( std::vector< T, A > const & a, vectorL< L, T, A > const & b );


/// @brief std::vector > vectorL
template< platform::SSize L, typename T, typename A >
bool
operator >( std::vector< T, A > const & a, vectorL< L, T, A > const & b );


/// @brief swap( vectorL, vectorL )
template< platform::SSize L, typename T, typename A >
void
swap( vectorL< L, T, A > & a, vectorL< L, T, A > & b );


/// @brief swap( vectorL, std::vector )
template< platform::SSize L, typename T, typename A >
void
swap( vectorL< L, T, A > & a, std::vector< T, A > & b );


/// @brief swap( std::vector, vectorL )
template< platform::SSize L, typename T, typename A >
void
swap( std::vector< T, A > & a, vectorL< L, T, A > & b );


} // namespace utility


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or utility::swap.  The legal alternative would be
// to add specializations of swap for each anticipated vectorL instantiation.


namespace std {


/// @brief swap( vectorL, vectorL )
template< platform::SSize L, typename T, typename A >
inline
void
swap( utility::vectorL< L, T, A > & a, utility::vectorL< L, T, A > & b )
{
	a.swap( b );
}


/// @brief swap( vectorL, std::vector )
template< platform::SSize L, typename T, typename A >
inline
void
swap( utility::vectorL< L, T, A > & a, std::vector< T, A > & b )
{
	a.swap( b );
}


/// @brief swap( std::vector, vectorL )
template< platform::SSize L, typename T, typename A >
inline
void
swap( std::vector< T, A > & a, utility::vectorL< L, T, A > & b )
{
	b.swap( a );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


// bool specialization
#include <utility/vectorL_bool.hh>


#endif // INCLUDED_utility_vectorL_HH
