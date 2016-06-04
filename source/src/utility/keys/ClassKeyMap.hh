// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/keys/ClassKeyMap.hh
/// @brief  Keyed-access map with key subset map for each client class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Key can be any type that is convertible to the index map's index type
///  @li The Key type should not be the same as Index or you'll get an operator[] ambiguity
///  @li If a utility Key subtype is used it must declare the ClassKeyMap as a friend
///  @li Keys are added to map by assign(), operator(), add(), and activate()
///  @li Keys can be added to map out of order
///  @li Client (C) parameter => Distinct type per client
///  @li Static index map => Client vectors share key set
///  @li Can create elements and then assign keys or vice versa: At any time you
///      Can have a vector with more or fewer elements than there are active keys


#ifndef INCLUDED_utility_keys_ClassKeyMap_hh
#define INCLUDED_utility_keys_ClassKeyMap_hh


// Unit headers
#include <utility/keys/ClassKeyMap.fwd.hh>

// Project headers
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>
#include <utility>


namespace utility {
namespace keys {


/// @brief Keyed-access map with key subset map for each client class
template< typename K, typename T, typename C >
class ClassKeyMap
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
	typedef  C  client_type;

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
	typedef  C  Client;


public: // Creation


	/// @brief Default constructor
	inline
	ClassKeyMap()
	{}


	/// @brief Copy constructor
	inline
	ClassKeyMap( ClassKeyMap const & a ) :
		v_( a.v_ )
	{}


	/// @brief Iterator range constructor
	template< typename InputIterator >
	inline
	ClassKeyMap(
		InputIterator const beg,
		InputIterator const end
	) :
		v_( beg, end )
	{
		for ( Index i = 1, e = v_.size(); i != e; ++i ) {
			add_key( v_[ i ].first );
		}
	}


	/// @brief Destructor
	inline
	~ClassKeyMap()
	{}


public: // Assignment


	/// @brief Copy assignment
	inline
	ClassKeyMap &
	operator =( ClassKeyMap const & a )
	{
		if ( this != &a ) {
			v_ = a.v_;
		}
		return *this;
	}


	/// @brief Uniform value assignment to current elements
	inline
	ClassKeyMap &
	operator =( Mapped const & mapped )
	{
		for ( Index i = 1, e = v_.size(); i <= e; ++i ) {
			v_[ i ].second = mapped;
		}
		return *this;
	}


	/// @brief Assign a value to an element with a key
	/// @note Adds the key to the map if not present
	/// @note Expands the vector if necessary
	inline
	ClassKeyMap &
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
	ClassKeyMap &
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
	ClassKeyMap &
	add( Key const & key )
	{
		add_key( key );
		return *this;
	}


	/// @brief Insert an element
	/// @note Adds the key to the map if not present
	/// @note Expands the vector if necessary
	inline
	ClassKeyMap &
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


	/// @brief Shrink the vector to remove unused capacity
	inline
	void
	shrink()
	{
		v_.shrink();
	}


	/// @brief swap( ClassKeyMap )
	inline
	void
	swap( ClassKeyMap & a )
	{
		v_.swap( a.v_ );
	}


	/// @brief swap( ClassKeyMap, ClassKeyMap )
	template< typename UK, typename UT, typename UC >
	friend
	void
	swap( ClassKeyMap< UK, UT, UC > & a, ClassKeyMap< UK, UT, UC > & b );

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


	/// @brief ClassKeyMap( key )
	/// @note Activates the key if inactive
	/// @note Expands the vector if necessary
	inline
	MappedReference
	operator ()( Key const & key )
	{
		return v_[ add_key( key ) ].second;
	}


	/// @brief ClassKeyMap[ key ] const
	inline
	MappedConstReference
	operator []( Key const & key ) const
	{
	debug_assert( active( key ) );
		return v_[ m()[ key ] ].second;
	}


	/// @brief ClassKeyMap[ key ]
	inline
	MappedReference
	operator []( Key const & key )
	{
	debug_assert( active( key ) );
		return v_[ m()[ key ] ].second;
	}


	/// @brief ClassKeyMap[ index ] const
	inline
	MappedConstReference
	operator []( Index const & i ) const
	{
		return v_[ i ].second;
	}


	/// @brief ClassKeyMap[ index ]
	inline
	MappedReference
	operator []( Index const & i )
	{
		return v_[ i ].second;
	}


	/// @brief ClassKeyMap( index ) const
	inline
	ConstReference
	operator ()( Index const & i ) const
	{
		return v_[ i ];
	}


	/// @brief ClassKeyMap( index )
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


	/// @brief ClassKeyMap == ClassKeyMap
	friend
	inline
	bool
	operator ==( ClassKeyMap const & a, ClassKeyMap const & b )
	{
		return ( a.v_ == b.v_ );
	}


	/// @brief ClassKeyMap != ClassKeyMap
	friend
	inline
	bool
	operator !=( ClassKeyMap const & a, ClassKeyMap const & b )
	{
		return ( a.v_ != b.v_ );
	}


	/// @brief ClassKeyMap < ClassKeyMap
	friend
	inline
	bool
	operator <( ClassKeyMap const & a, ClassKeyMap const & b )
	{
		return ( a.v_ < b.v_ );
	}


	/// @brief ClassKeyMap > ClassKeyMap
	friend
	inline
	bool
	operator >( ClassKeyMap const & a, ClassKeyMap const & b )
	{
		return ( a.v_ > b.v_ );
	}


	/// @brief ClassKeyMap <= ClassKeyMap
	friend
	inline
	bool
	operator <=( ClassKeyMap const & a, ClassKeyMap const & b )
	{
		return ( a.v_ <= b.v_ );
	}


	/// @brief ClassKeyMap >= ClassKeyMap
	friend
	inline
	bool
	operator >=( ClassKeyMap const & a, ClassKeyMap const & b )
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


}; // ClassKeyMap


// Friend function namespace declarations


/// @brief swap( ClassKeyMap, ClassKeyMap )
template< typename K, typename T, typename C >
void
swap( ClassKeyMap< K, T, C > & a, ClassKeyMap< K, T, C > & b )
{
	a.v_.swap( b.v_ );
}


/// @brief ClassKeyMap == ClassKeyMap
template< typename K, typename T, typename C >
bool
operator ==( ClassKeyMap< K, T, C > const & a, ClassKeyMap< K, T, C > const & b );


/// @brief ClassKeyMap != ClassKeyMap
template< typename K, typename T, typename C >
bool
operator !=( ClassKeyMap< K, T, C > const & a, ClassKeyMap< K, T, C > const & b );


/// @brief ClassKeyMap < ClassKeyMap
template< typename K, typename T, typename C >
bool
operator <( ClassKeyMap< K, T, C > const & a, ClassKeyMap< K, T, C > const & b );


/// @brief ClassKeyMap > ClassKeyMap
template< typename K, typename T, typename C >
bool
operator >( ClassKeyMap< K, T, C > const & a, ClassKeyMap< K, T, C > const & b );


/// @brief ClassKeyMap <= ClassKeyMap
template< typename K, typename T, typename C >
bool
operator <=( ClassKeyMap< K, T, C > const & a, ClassKeyMap< K, T, C > const & b );


/// @brief ClassKeyMap >= ClassKeyMap
template< typename K, typename T, typename C >
bool
operator >=( ClassKeyMap< K, T, C > const & a, ClassKeyMap< K, T, C > const & b );


} // namespace keys
} // namespace utility


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add template functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or utility::swap.  The legal alternative would be
// to add specializations of swap for each anticipated ClassKeyMap instantiation.


namespace std {


/// @brief swap( ClassKeyMap, ClassKeyMap )
template< typename K, typename T, typename C >
inline
void
swap( utility::keys::ClassKeyMap< K, T, C > & a, utility::keys::ClassKeyMap< K, T, C > & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_utility_keys_ClassKeyMap_HH
