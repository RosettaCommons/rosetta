// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/keys/SmallKeyVector.hh
/// @brief  Keyed-access vector with key subset map
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Key can be any type that is convertible to the index map's index type
///  @li The Key type should not be the same as Index or you'll get an operator[] ambiguity
///  @li If a utility Key subtype is used it must declare the SmallKeyVector as a friend
///  @li Keys are added to map by assign(), operator(), and add()
///  @li Keys can be added to map out of order
///  @li Index map is specific to the SmallKeyVector so this container is intended for use
///      where the index range of the keys is either small or dense
///  @li Can create elements and then assign keys or vice versa: At any time you
///      Can have a vector with more or fewer elements than there are active keys


#ifndef INCLUDED_utility_keys_SmallKeyVector_hh
#define INCLUDED_utility_keys_SmallKeyVector_hh


// Unit headers
#include <utility/keys/SmallKeyVector.fwd.hh>

// Project headers
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>
#include <utility/assert.hh>


namespace utility {
namespace keys {


/// @brief Keyed-access vector with key subset map
template< typename K, typename T >
class SmallKeyVector
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


public: // Creation


	/// @brief Default constructor
	inline
	SmallKeyVector() :
		u_( 0 )
	{}


	/// @brief Copy constructor
	inline
	SmallKeyVector( SmallKeyVector const & a ) :
		v_( a.v_ ),
		m_( a.m_ ),
		u_( a.u_ )
	{}


	/// @brief Size constructor
	inline
	explicit
	SmallKeyVector( Size const num ) :
		v_( num ),
		u_( 0 )
	{}


	/// @brief Uniform value constructor
	inline
	SmallKeyVector(
		Size const num,
		Value const & value
	) :
		v_( num, value ),
		u_( 0 )
	{}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	SmallKeyVector(
		InputIterator const beg,
		InputIterator const end
	) :
		v_( beg, end ),
		u_( 0 )
	{}


	/// @brief Destructor
	inline
	~SmallKeyVector()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	SmallKeyVector &
	operator =( SmallKeyVector const & a )
	{
		if ( this != &a ) {
			v_ = a.v_;
			m_ = a.m_;
			u_ = a.u_;
		}
		return *this;
	}


	/// @brief Uniform value assignment to current elements
	inline
	SmallKeyVector &
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
	SmallKeyVector &
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
	SmallKeyVector &
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


	/// @brief Shrink the vectors to remove unused capacity
	inline
	void
	shrink()
	{
		v_.shrink();
		m_.shrink();
	}


	/// @brief swap( SmallKeyVector )
	inline
	void
	swap( SmallKeyVector & a )
	{
		v_.swap( a.v_ );
		m_.swap( a.m_ );
		std::swap( u_, a.u_ );
	}


	/// @brief swap( SmallKeyVector, SmallKeyVector )
	friend
	inline
	void
	swap( SmallKeyVector & a, SmallKeyVector & b )
	{
		a.v_.swap( b.v_ );
		a.m_.swap( b.m_ );
		std::swap( a.u_, b.u_ );
	}


	/// @brief Clear the vector
	inline
	void
	clear()
	{
		v_.clear();
		m_.clear();
		u_ = 0;
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
		return active( key );
	}


	/// @brief Is a key active?
	inline
	bool
	active( Key const & key ) const
	{
		return ( ( m_.has( key ) ) && ( m_[ key ] != 0 ) );
	}


	/// @brief Is a key inactive?
	inline
	bool
	inactive( Key const & key ) const
	{
		return ( ( ! m_.has( key ) ) || ( m_[ key ] == 0 ) );
	}


	/// @brief Index of a key
	inline
	Index const &
	index( Key const & key )
	{
		debug_assert( active( key ) );
		return m_[ key ];
	}


	/// @brief Iterator to element with a key
	inline
	ConstIterator
	find( Key const & key ) const
	{
		return ( active( key ) ? v_.begin() + m_[ key ] - 1 : v_.end() );
	}


	/// @brief Iterator to element with a key
	inline
	Iterator
	find( Key const & key )
	{
		return ( active( key ) ? v_.begin() + m_[ key ] - 1 : v_.end() );
	}


public: // Indexers


	/// @brief SmallKeyVector( key )
	/// @note Activates the key if inactive
	/// @note Expands the vector if necessary
	inline
	Reference
	operator ()( Key const & key )
	{
		return v_[ add_key( key ) ];
	}


	/// @brief SmallKeyVector[ key ] const
	inline
	ConstReference
	operator []( Key const & key ) const
	{
		debug_assert( active( key ) );
		return v_[ m_[ key ] ];
	}


	/// @brief SmallKeyVector[ key ]
	inline
	Reference
	operator []( Key const & key )
	{
		debug_assert( active( key ) );
		return v_[ m_[ key ] ];
	}


	/// @brief Element at index key: Bounds checked
	inline
	ConstReference
	at( Key const & key ) const
	{
		debug_assert( active( key ) );
		return v_.at( m_.at( key ) );
	}


	/// @brief Element at index key: Bounds checked
	inline
	Reference
	at( Key const & key )
	{
		debug_assert( active( key ) );
		return v_.at( m_.at( key ) );
	}


	/// @brief SmallKeyVector[ index ] const
	inline
	ConstReference
	operator []( Index const & i ) const
	{
		return v_[ i ];
	}


	/// @brief SmallKeyVector[ index ]
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


	/// @brief SmallKeyVector == SmallKeyVector
	friend
	inline
	bool
	operator ==( SmallKeyVector const & a, SmallKeyVector const & b )
	{
		return ( ( a.v_ == b.v_ ) && ( a.m_ == b.m_ ) );
	}


	/// @brief SmallKeyVector != SmallKeyVector
	friend
	inline
	bool
	operator !=( SmallKeyVector const & a, SmallKeyVector const & b )
	{
		return !( a == b );
	}


private: // Methods


	/// @brief Add an element with a key if not present and return its index: Activate key if inactive
	inline
	Index const &
	add_key( Key const & key )
	{
		debug_assert( key > 0 );
		if ( ! m_.has( key ) ) { // Extend index map and activate key
			m_.resize( key, Index( 0 ) );
			m_[ key ] = ++u_;
		} else if ( m_[ key ] == 0 ) { // Activate key
			m_[ key ] = ++u_;
		}
		Index const & i( m_[ key ] );
		if ( i > v_.size() ) v_.resize( i );
		return i;
	}


private: // Fields


	/// @brief Vector of values indexed by a subset of the possible keys
	Vector v_;

	/// @brief Index map from keys into v_: Zero => inactive key
	IndexMap m_;

	/// @brief Upper active index of active keys
	Index u_;


}; // SmallKeyVector


// Friend function namespace declarations


/// @brief swap( SmallKeyVector, SmallKeyVector )
template< typename K, typename T >
void
swap( SmallKeyVector< K, T > & a, SmallKeyVector< K, T > & b );


/// @brief SmallKeyVector == SmallKeyVector
template< typename K, typename T >
bool
operator ==( SmallKeyVector< K, T > const & a, SmallKeyVector< K, T > const & b );


/// @brief SmallKeyVector != SmallKeyVector
template< typename K, typename T >
bool
operator !=( SmallKeyVector< K, T > const & a, SmallKeyVector< K, T > const & b );


} // namespace keys
} // namespace utility


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or utility::swap.  The legal alternative would be
// to add specializations of swap for each anticipated SmallKeyVector instantiation.


namespace std {


/// @brief swap( SmallKeyVector, SmallKeyVector )
template< typename K, typename T >
inline
void
swap( utility::keys::SmallKeyVector< K, T > & a, utility::keys::SmallKeyVector< K, T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_utility_keys_SmallKeyVector_HH
