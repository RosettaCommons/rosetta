// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/SmallKeyMap.hh
/// @brief  Keyed-access map with key subset map
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Key can be any type that is convertible to the index map's index type
///  @li The Key type should not be the same as Index or you'll get an operator[] ambiguity
///  @li If a utility Key subtype is used it must declare the SmallKeyMap as a friend
///  @li Keys are added to map by assign(), operator(), and add()
///  @li Keys can be added to map out of order
///  @li Index map is specific to the SmallKeyMap so this container is intended for use
///      where the index range of the keys is either small or dense


#ifndef INCLUDED_utility_keys_SmallKeyMap_hh
#define INCLUDED_utility_keys_SmallKeyMap_hh


// Unit headers
#include <utility/keys/SmallKeyMap.fwd.hh>

// Project headers
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>
#include <utility/assert.hh>
#include <utility>


namespace utility {
namespace keys {


/// @brief Keyed-access map with key subset map
template< typename K, typename T >
class SmallKeyMap
{


private: // Types


	typedef  vector1< std::pair< K, T > >  Vector;
	typedef  vector1< typename Vector::Index >  IndexMap;
	typedef  typename IndexMap::Size  IndexMapSize;
	typedef  typename IndexMap::Index  IndexMapIndex;


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  T  mapped_type;
	typedef  T &  mapped_reference;
	typedef  T const &  mapped_const_reference;
	typedef  T *  mapped_pointer;
	typedef  T const *  mapped_const_pointer;
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
	typedef  T  Mapped;
	typedef  T &  MappedReference;
	typedef  T const &  MappedConstReference;
	typedef  T *  MappedPointer;
	typedef  T const *  MappedConstPointer;
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
	SmallKeyMap() :
		u_( 0 )
	{}


	/// @brief Copy constructor
	inline
	SmallKeyMap( SmallKeyMap const & a ) :
		v_( a.v_ ),
		m_( a.m_ ),
		u_( a.u_ )
	{}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	SmallKeyMap(
		InputIterator const beg,
		InputIterator const end
	) :
		v_( beg, end ),
		u_( 0 )
	{
		for ( Index i = 1, e = v_.size(); i != e; ++i ) {
			add_key( v_[ i ].first );
		}
	}


	/// @brief Destructor
	inline
	~SmallKeyMap()
	= default;


public: // Assignment


	/// @brief Copy assignment
	inline
	SmallKeyMap &
	operator =( SmallKeyMap const & a )
	{
		if ( this != &a ) {
			v_ = a.v_;
			m_ = a.m_;
			u_ = a.u_;
		}
		return *this;
	}


	/// @brief Uniform mapped value assignment to current elements
	inline
	SmallKeyMap &
	operator =( Mapped const & mapped )
	{
		for ( Index i = 1, e = v_.size(); i <= e; ++i ) {
			v_[ i ].second = mapped;
		}
		return *this;
	}


	/// @brief Assign a mapped value to an element with a key
	/// @note Adds the key to the map if not present
	/// @note Expands the vector if necessary
	inline
	SmallKeyMap &
	assign(
		Key const & key,
		Mapped const & mapped
	)
	{
		v_[ add_key( key ) ] = Value( key, mapped );
		return *this;
	}


	/// @brief Assign a value to an element
	/// @note Adds the key to the map if not present
	/// @note Expands the vector if necessary
	inline
	SmallKeyMap &
	assign( Value const & value )
	{
		v_[ add_key( value.first ) ] = value;
		return *this;
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
		clear();
		v_.assign( beg, end );
		for ( Index i = 1, e = v_.size(); i != e; ++i ) {
			add_key( v_[ i ].first );
		}
	}


public: // Methods


	/// @brief Add an element with a key if not present: Activate key if inactive
	inline
	SmallKeyMap &
	add( Key const & key )
	{
		add_key( key );
		return *this;
	}


	/// @brief Insert an element
	/// @note Adds the key to the map if not present
	/// @note Expands the vector if necessary
	inline
	SmallKeyMap &
	insert( Value const & value )
	{
		v_[ add_key( value.first ) ] = value;
		return *this;
	}


	/// @brief Insert elements from iterator range [beg,end)
	template< typename InputIterator >
	inline
	void
	insert(
		InputIterator const beg,
		InputIterator const end
	)
	{
		for ( InputIterator i = beg; i != end; ++i ) {
			insert( *i );
		}
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


	/// @brief swap( SmallKeyMap )
	inline
	void
	swap( SmallKeyMap & a )
	{
		v_.swap( a.v_ );
		m_.swap( a.m_ );
		std::swap( u_, a.u_ );
	}


	/// @brief swap( SmallKeyMap, SmallKeyMap )
	friend
	inline
	void
	swap( SmallKeyMap & a, SmallKeyMap & b )
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


	/// @brief SmallKeyMap( key )
	/// @note Activates the key if inactive
	/// @note Expands the vector if necessary
	inline
	MappedReference
	operator ()( Key const & key )
	{
		return v_[ add_key( key ) ].second;
	}


	/// @brief SmallKeyMap[ key ] const
	inline
	MappedConstReference
	operator []( Key const & key ) const
	{
	debug_assert( active( key ) );
		return v_[ m_[ key ] ].second;
	}


	/// @brief SmallKeyMap[ key ]
	inline
	MappedReference
	operator []( Key const & key )
	{
	debug_assert( active( key ) );
		return v_[ m_[ key ] ].second;
	}


	/// @brief SmallKeyMap[ index ] const
	inline
	MappedConstReference
	operator []( Index const & i ) const
	{
		return v_[ i ].second;
	}


	/// @brief SmallKeyMap[ index ]
	inline
	MappedReference
	operator []( Index const & i )
	{
		return v_[ i ].second;
	}


	/// @brief SmallKeyMap( index ) const
	inline
	ConstReference
	operator ()( Index const & i ) const
	{
		return v_[ i ];
	}


	/// @brief SmallKeyMap( index )
	inline
	Reference
	operator ()( Index const & i )
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


	/// @brief SmallKeyMap == SmallKeyMap
	friend
	inline
	bool
	operator ==( SmallKeyMap const & a, SmallKeyMap const & b )
	{
		return ( ( a.v_ == b.v_ ) && ( a.m_ == b.m_ ) );
	}


	/// @brief SmallKeyMap != SmallKeyMap
	friend
	inline
	bool
	operator !=( SmallKeyMap const & a, SmallKeyMap const & b )
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


}; // SmallKeyMap


// Friend function namespace declarations


/// @brief swap( SmallKeyMap, SmallKeyMap )
template< typename K, typename T >
void
swap( SmallKeyMap< K, T > & a, SmallKeyMap< K, T > & b );


/// @brief SmallKeyMap == SmallKeyMap
template< typename K, typename T >
bool
operator ==( SmallKeyMap< K, T > const & a, SmallKeyMap< K, T > const & b );


/// @brief SmallKeyMap != SmallKeyMap
template< typename K, typename T >
bool
operator !=( SmallKeyMap< K, T > const & a, SmallKeyMap< K, T > const & b );


} // namespace keys
} // namespace utility


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or utility::swap.  The legal alternative would be
// to add specializations of swap for each anticipated SmallKeyMap instantiation.


namespace std {


/// @brief swap( SmallKeyMap, SmallKeyMap )
template< typename K, typename T >
inline
void
swap( utility::keys::SmallKeyMap< K, T > & a, utility::keys::SmallKeyMap< K, T > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_utility_keys_SmallKeyMap_HH
