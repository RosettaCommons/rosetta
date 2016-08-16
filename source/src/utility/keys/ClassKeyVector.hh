// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/ClassKeyVector.hh
/// @brief  Keyed-access vector with key subset map for each client class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Key can be any type that is convertible to the index map's index type
///  @li The Key type should not be the same as Index or you'll get an operator[] ambiguity
///  @li If a utility Key subtype is used it must declare the ClassKeyVector as a friend
///  @li Keys are added to map by assign(), operator(), add(), and activate()
///  @li Keys can be added to map out of order
///  @li Client (C) parameter => Distinct type per client
///  @li Static index map => Client vectors share key set
///  @li Can create elements and then assign keys or vice versa: At any time you
///      Can have a vector with more or fewer elements than there are active keys


#ifndef INCLUDED_utility_keys_ClassKeyVector_hh
#define INCLUDED_utility_keys_ClassKeyVector_hh


// Unit headers
#include <utility/keys/ClassKeyVector.fwd.hh>

// Project headers
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>
#include <iterator>
#include <algorithm> //Required by GCC 4.3.2

namespace utility {
namespace keys {


/// @brief Keyed-access vector with key subset map for each client class
template< typename K, typename T, typename C >
class ClassKeyVector
{


private: // Types


	typedef  vector1< T >  Vector;
	typedef  vector1< typename Vector::Index >  IndexMap;
	typedef  typename IndexMap::Size  IndexMapSize;
	typedef  typename IndexMap::Index  IndexMapIndex;


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  typename Vector::value_type  value_type;
	typedef  typename Vector::reference  reference;
	typedef  typename Vector::const_reference  const_reference;
	typedef  typename Vector::pointer  pointer;
	typedef  typename Vector::const_pointer  const_pointer;
	typedef  typename Vector::iterator  iterator;
	typedef  typename Vector::const_iterator  const_iterator;
	typedef  typename Vector::reverse_iterator  reverse_iterator;
	typedef  typename Vector::const_reverse_iterator  const_reverse_iterator;
	typedef  typename Vector::size_type  size_type;
	typedef  typename Vector::index_type  index_type;
	typedef  typename Vector::difference_type  difference_type;
	typedef  typename Vector::allocator_type  allocator_type;
	typedef  C  client_type;

	// Project style
	typedef  K  Key;
	typedef  typename Vector::Value  Value;
	typedef  typename Vector::Reference  Reference;
	typedef  typename Vector::ConstReference  ConstReference;
	typedef  typename Vector::Pointer  Pointer;
	typedef  typename Vector::ConstPointer  ConstPointer;
	typedef  typename Vector::Iterator  Iterator;
	typedef  typename Vector::ConstIterator  ConstIterator;
	typedef  typename Vector::ReverseIterator  ReverseIterator;
	typedef  typename Vector::ConstReverseIterator  ConstReverseIterator;
	typedef  typename Vector::Size  Size;
	typedef  typename Vector::Index  Index;
	typedef  typename Vector::Difference  Difference;
	typedef  typename Vector::Allocator  Allocator;
	typedef  C  Client;


public: // Creation


	/// @brief Default constructor
	inline
	ClassKeyVector()
	{}


	/// @brief Copy constructor
	inline
	ClassKeyVector( ClassKeyVector const & a ) :
		v_( a.v_ )
	{}


	/// @brief Size constructor
	inline
	explicit
	ClassKeyVector( Size const num ) :
		v_( num )
	{}


	/// @brief Uniform value constructor
	inline
	ClassKeyVector(
		Size const num,
		Value const & value
	) :
		v_( num, value )
	{}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	ClassKeyVector(
		InputIterator const beg,
		InputIterator const end
	) :
		v_( beg, end )
	{}


	/// @brief Destructor
	inline
	~ClassKeyVector()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	ClassKeyVector &
	operator =( ClassKeyVector const & a )
	{
		if ( this != &a ) {
			v_ = a.v_;
		}
		return *this;
	}


	/// @brief Uniform value assignment to current elements
	inline
	ClassKeyVector &
	operator =( Value const & value )
	{
		for ( Index i = 1, e = v_.size(); i <= e; ++i ) {
			v_[ i ] = value;
		}
		return *this;
	}


	/// @brief Assign a value to an element with a key
	/// @note Adds the key to the map if not present
	/// @note Expands the vector if necessary
	inline
	ClassKeyVector &
	assign(
		Key const & key,
		Value const & value
	)
	{
		v_[ add_key( key ) ] = value;
		return *this;
	}


	/// @brief Uniform value assignment
	inline
	void
	assign(
		Size const num,
		Value const & value
	)
	{
		v_.assign( num, value );
	}


	/// @brief Iterator assignment
	template< typename InputIterator >
	inline
	void
	assign(
		InputIterator const beg,
		InputIterator const end
	)
	{
		v_.assign( beg, end );
	}


public: // Methods


	/// @brief Add an element with a key if not present: Activate key if inactive
	inline
	ClassKeyVector &
	add( Key const & key )
	{
		add_key( key );
		return *this;
	}


	/// @brief Insert an element at an iterator position
	inline
	Iterator
	insert(
		Iterator const pos,
		Value const & value
	)
	{
		return v_.insert( pos, value );
	}


	/// @brief Insert num copies of an element at an iterator position
	inline
	void
	insert(
		Iterator const pos,
		Size const num,
		Value const & value
	)
	{
		v_.insert( pos, num, value );
	}


	/// @brief Insert elements from iterator range [beg,end) at an iterator position
	template< typename InputIterator >
	inline
	void
	insert(
		Iterator const pos,
		InputIterator const beg,
		InputIterator const end
	)
	{
		v_.insert( pos, beg, end );
	}


	/// @brief Appends an element
	inline
	void
	push_back( Value const & value )
	{
		v_.push_back( value );
	}


	/// @brief Erase an element at an iterator position
	inline
	Iterator
	erase( Iterator const pos )
	{
		return v_.erase( pos );
	}


	/// @brief Erase elements in the iterator range [beg,end)
	inline
	Iterator
	erase(
		Iterator const beg,
		Iterator const end
	)
	{
		return v_.erase( beg, end );
	}


	/// @brief Removes the last element
	inline
	void
	pop_back()
	{
	debug_assert( ! v_.empty() );
		v_.pop_back();
	}


	/// @brief Resize: Default construct new elements
	inline
	void
	resize( Size const num )
	{
		v_.resize( num );
	}


	/// @brief Resize: Assign given value to new elements
	inline
	void
	resize(
		Size const num,
		Value const & value
	)
	{
		v_.resize( num, value );
	}


	/// @brief Reserve space for a given number of elements
	inline
	void
	reserve( Size const num )
	{
		v_.reserve( num );
	}


	/// @brief Shrink the vector to remove unused capacity
	inline
	void
	shrink()
	{
		v_.shrink();
	}


	/// @brief Sort the vector with a given predicate
	template< typename Compare >
	inline
	void
	sort( Compare compare )
	{
		Vector const s( v_ ); // Save original vector
		std::sort( v_.begin(), v_.end(), compare );
		if ( s != v_ ) { // Adjust index map for sort
			IndexMap & imap( m() );
			for ( typename IndexMap::iterator k = imap.begin(), ke = imap.end(); k != ke; ++k ) {
				Index & i( *k ); // Old index
				if ( i != 0 ) { // Active key
					i = std::distance( v_.begin(), std::lower_bound( v_.begin(), v_.end(), s[ i ], compare ) ) + 1; // New index
				}
			}
		}
	}


	/// @brief swap( ClassKeyVector )
	inline
	void
	swap( ClassKeyVector & a )
	{
		v_.swap( a.v_ );
	}


	/// @brief swap( ClassKeyVector, ClassKeyVector )
	template< typename UK, typename UT, typename UC >
	friend
	void
	swap( ClassKeyVector< UK, UT, UC > & a, ClassKeyVector< UK, UT, UC > & b );


	/// @brief Clear the vector
	inline
	void
	clear()
	{
		v_.clear();
	}


public: // Properties


	/// @brief Size
	inline
	Size
	size() const
	{
		return v_.size();
	}


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return v_.empty();
	}


	/// @brief Max size
	inline
	Size
	max_size() const
	{
		return v_.max_size();
	}


	/// @brief Capacity
	inline
	Size
	capacity() const
	{
		return v_.capacity();
	}


	/// @brief Front element
	inline
	ConstReference
	front() const
	{
	debug_assert( ! v_.empty() );
		return v_.front();
	}


	/// @brief Front element
	inline
	Reference
	front()
	{
	debug_assert( ! v_.empty() );
		return v_.front();
	}


	/// @brief Back element
	inline
	ConstReference
	back() const
	{
	debug_assert( ! v_.empty() );
		return v_.back();
	}


	/// @brief Back element
	inline
	Reference
	back()
	{
	debug_assert( ! v_.empty() );
		return v_.back();
	}


	/// @brief Is an element with a key present?
	inline
	bool
	has( Key const & key ) const
	{
		return ( ( key > 0 ) && ( m().has( key ) ) && ( v_.has( m()[ key ] ) ) );
	}


	/// @brief Index of a key
	inline
	static
	Index const &
	index( Key const & key )
	{
	debug_assert( active( key ) );
		return m()[ key ];
	}


	/// @brief Iterator to element with a key
	inline
	ConstIterator
	find( Key const & key ) const
	{
		return ( active( key ) ? v_.begin() + m()[ key ] - 1 : v_.end() );
	}


	/// @brief Iterator to element with a key
	inline
	Iterator
	find( Key const & key )
	{
		return ( active( key ) ? v_.begin() + m()[ key ] - 1 : v_.end() );
	}


public: // Indexers


	/// @brief ClassKeyVector( key )
	/// @note Activates the key if inactive
	/// @note Expands the vector if necessary
	inline
	Reference
	operator ()( Key const & key )
	{
		return v_[ add_key( key ) ];
	}


	/// @brief ClassKeyVector[ key ] const
	inline
	ConstReference
	operator []( Key const & key ) const
	{
	debug_assert( active( key ) );
		return v_[ m()[ key ] ];
	}


	/// @brief ClassKeyVector[ key ]
	inline
	Reference
	operator []( Key const & key )
	{
	debug_assert( active( key ) );
		return v_[ m()[ key ] ];
	}


	/// @brief Element at index key: Bounds checked
	inline
	ConstReference
	at( Key const & key ) const
	{
	debug_assert( active( key ) );
		return v_.at( m().at( key ) );
	}


	/// @brief Element at index key: Bounds checked
	inline
	Reference
	at( Key const & key )
	{
	debug_assert( active( key ) );
		return v_.at( m().at( key ) );
	}


	/// @brief ClassKeyVector[ index ] const
	inline
	ConstReference
	operator []( Index const & i ) const
	{
		return v_[ i ];
	}


	/// @brief ClassKeyVector[ index ]
	inline
	Reference
	operator []( Index const & i )
	{
		return v_[ i ];
	}


public: // Iterators


	/// @brief Begin iterator
	inline
	ConstIterator
	begin() const
	{
		return v_.begin();
	}


	/// @brief Begin iterator
	inline
	Iterator
	begin()
	{
		return v_.begin();
	}


	/// @brief End iterator
	inline
	ConstIterator
	end() const
	{
		return v_.end();
	}


	/// @brief End iterator
	inline
	Iterator
	end()
	{
		return v_.end();
	}


	/// @brief Begin reverse iterator
	inline
	ConstReverseIterator
	rbegin() const
	{
		return v_.rbegin();
	}


	/// @brief Begin reverse iterator
	inline
	ReverseIterator
	rbegin()
	{
		return v_.rbegin();
	}


	/// @brief End reverse iterator
	inline
	ConstReverseIterator
	rend() const
	{
		return v_.rend();
	}


	/// @brief End reverse iterator
	inline
	ReverseIterator
	rend()
	{
		return v_.rend();
	}


public: // Comparison


	/// @brief ClassKeyVector == ClassKeyVector
	friend
	inline
	bool
	operator ==( ClassKeyVector const & a, ClassKeyVector const & b )
	{
		return ( a.v_ == b.v_ );
	}


	/// @brief ClassKeyVector != ClassKeyVector
	friend
	inline
	bool
	operator !=( ClassKeyVector const & a, ClassKeyVector const & b )
	{
		return ( a.v_ != b.v_ );
	}


	/// @brief ClassKeyVector < ClassKeyVector
	friend
	inline
	bool
	operator <( ClassKeyVector const & a, ClassKeyVector const & b )
	{
		return ( a.v_ < b.v_ );
	}


	/// @brief ClassKeyVector > ClassKeyVector
	friend
	inline
	bool
	operator >( ClassKeyVector const & a, ClassKeyVector const & b )
	{
		return ( a.v_ > b.v_ );
	}


	/// @brief ClassKeyVector <= ClassKeyVector
	friend
	inline
	bool
	operator <=( ClassKeyVector const & a, ClassKeyVector const & b )
	{
		return ( a.v_ <= b.v_ );
	}


	/// @brief ClassKeyVector >= ClassKeyVector
	friend
	inline
	bool
	operator >=( ClassKeyVector const & a, ClassKeyVector const & b )
	{
		return ( a.v_ >= b.v_ );
	}


public: // Static functions


	/// @brief Is a key active?
	inline
	static
	bool
	active( Key const & key )
	{
		return ( ( m().has( key ) ) && ( m()[ key ] != 0 ) );
	}


	/// @brief Is a key inactive?
	inline
	static
	bool
	inactive( Key const & key )
	{
		return ( ( ! m().has( key ) ) || ( m()[ key ] == 0 ) );
	}


	/// @brief Activate a key if inactive
	inline
	static
	void
	activate( Key const & key )
	{
	debug_assert( key > 0 );
		IndexMap & imap( m() );
		if ( ! imap.has( key ) ) { // Extend index map and activate key
			imap.resize( key, Index( 0 ) );
			imap[ key ] = ++u();
		} else if ( imap[ key ] == 0 ) { // Activate key
			imap[ key ] = ++u();
		}
	}


private: // Methods


	/// @brief Add an element with a key if not present and return its index: Activate key if inactive
	inline
	Index const &
	add_key( Key const & key )
	{
	debug_assert( key > 0 );
		Index const & i( activated_index( key ) );
		if ( i > v_.size() ) v_.resize( i );
		return i;
	}


private: // Static functions


	/// @brief Index map from keys into v_: Zero => inactive key
	inline
	static
	IndexMap &
	m()
	{
		static IndexMap m_; // Function local to support globals
		return m_;
	}


	/// @brief Upper active index of active keys
	inline
	static
	Index &
	u()
	{
		static Index u_( 0 ); // Function local to support globals
		return u_;
	}


	/// @brief Activate a key if inactive and return its index
	inline
	static
	Index &
	activated_index( Key const & key )
	{
		activate( key );
		return m()[ key ];
	}


	/// @brief Shrink the index map to remove unused capacity
	inline
	static
	void
	map_shrink()
	{
		m().shrink();
	}


private: // Fields


	/// @brief Vector of values indexed by a subset of the possible keys
	Vector v_;


}; // ClassKeyVector


// Friend function namespace declarations


/// @brief swap( ClassKeyVector, ClassKeyVector )
template< typename K, typename T, typename C >
void
swap( ClassKeyVector< K, T, C > & a, ClassKeyVector< K, T, C > & b )
{
	a.v_.swap( b.v_ );
}


/// @brief ClassKeyVector == ClassKeyVector
template< typename K, typename T, typename C >
bool
operator ==( ClassKeyVector< K, T, C > const & a, ClassKeyVector< K, T, C > const & b );


/// @brief ClassKeyVector != ClassKeyVector
template< typename K, typename T, typename C >
bool
operator !=( ClassKeyVector< K, T, C > const & a, ClassKeyVector< K, T, C > const & b );


/// @brief ClassKeyVector < ClassKeyVector
template< typename K, typename T, typename C >
bool
operator <( ClassKeyVector< K, T, C > const & a, ClassKeyVector< K, T, C > const & b );


/// @brief ClassKeyVector > ClassKeyVector
template< typename K, typename T, typename C >
bool
operator >( ClassKeyVector< K, T, C > const & a, ClassKeyVector< K, T, C > const & b );


/// @brief ClassKeyVector <= ClassKeyVector
template< typename K, typename T, typename C >
bool
operator <=( ClassKeyVector< K, T, C > const & a, ClassKeyVector< K, T, C > const & b );


/// @brief ClassKeyVector >= ClassKeyVector
template< typename K, typename T, typename C >
bool
operator >=( ClassKeyVector< K, T, C > const & a, ClassKeyVector< K, T, C > const & b );


} // namespace keys
} // namespace utility


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or utility::swap.  The legal alternative would be
// to add specializations of swap for each anticipated ClassKeyVector instantiation.


namespace std {


/// @brief swap( ClassKeyVector, ClassKeyVector )
template< typename K, typename T, typename C >
inline
void
swap( utility::keys::ClassKeyVector< K, T, C > & a, utility::keys::ClassKeyVector< K, T, C > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_utility_keys_ClassKeyVector_HH
