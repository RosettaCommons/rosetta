// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/KeyLookup.hh
/// @brief  Key lookup map and collection and functors
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Monostate singleton: map and collection are static
///  @li Key (K) parameter: The key type
///  @li Client (C) parameter: The client (user) of the KeyLookup
///  @li Key must provide id(), identifier(), and code() members returning std::string
///  @li KeyLookup doesn't use map::operator[] to avoid default constructing keys,
///      which would trigger an infinite cycle of add/insert calls
///  @li Holds copy of key in VariantKey to preserve dynamic key type
///  @li Holds just one copy of each key with a given value (determined by its operator <)
///      even if the identifiers vary
///  @li If same value key is added again with any different identifiers all the
///      identifiers refer to the original key in the lookup map
///  @li An identifier can only be used for a single key value: assert will catch violation
///  @li Keys are cloned by VariantKey so add them in their actual type constructor body
///  @li Functors for KeyLookup queries and key iterators are provided for key collections


#ifndef INCLUDED_utility_keys_KeyLookup_hh
#define INCLUDED_utility_keys_KeyLookup_hh


// Unit headers
#include <utility/keys/KeyLookup.fwd.hh>

// Package headers
#include <utility/keys/Key.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/VariantKey.hh>

// C++ headers
#include <utility/assert.hh>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <algorithm>


namespace utility {
namespace keys {


/// @brief Key lookup map and collection
template< typename K >
class KeyLookup
{


private: // Types


	typedef  VariantKey< K >  Key_;
	typedef  std::set< Key_ >  Keys;
	typedef  Key_ const *  KeyP;
	typedef  std::map< std::string, KeyP >  Map;


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  typename Map::mapped_type  mapped_type;
	typedef  typename Map::value_type  value_type;
	typedef  typename Map::reference  reference;
	typedef  typename Map::const_reference  const_reference;
	typedef  typename Map::size_type  size_type;
	typedef  typename Keys::const_iterator  const_iterator;

	// Project style
	typedef  K  Key;
	typedef  typename Map::mapped_type  Mapped;
	typedef  typename Map::value_type  Value;
	typedef  typename Map::reference  Reference;
	typedef  typename Map::const_reference  ConstReference;
	typedef  typename Map::size_type  Size;
	typedef  typename Keys::const_iterator  ConstIterator;


private: // Creation


	/// @brief Default constructor
	inline
	KeyLookup()
	= default;


	/// @brief Copy constructor
	inline
	KeyLookup( KeyLookup const & a ); // Undefined


	/// @brief Destructor
	inline
	~KeyLookup()
	= default;


private: // Assignment


	/// @brief Copy assignment
	inline
	KeyLookup &
	operator =( KeyLookup const & a ); // Undefined


public: // Methods


	/// @brief Add a key
	inline
	static
	void
	add( Key const & key )
	{
		// Put key in collection
		Key_ const & key_( *(keys().insert( key ).first) );

		// Put key in lookup map with all three identifiers
		insert( key.id(), key_ );
		insert( key.identifier(), key_ );
		insert( key.code(), key_ );
	}


	/// @brief Add a key with an extra (non-unique) identifier
	inline
	static
	void
	add( Key const & key, std::string const & identifier )
	{
		// Put key in collection
		Key_ const & key_( *(keys().insert( key ).first) );

		// Put key in lookup map with all four identifiers
		insert( key.id(), key_ );
		insert( key.identifier(), key_ );
		insert( key.code(), key_ );
		insert_nonunique( identifier, key_ );
	}


	/// @brief Key with an identifier: Generate key if not present
	template< typename KeyC > // Concrete Key type to create
	inline
	static
	Key const &
	gen( std::string const & id )
	{
		Map const & m_( m() );
		typename Map::const_iterator const i( m_.find( stripped_whitespace( id ) ) );
		if ( i != m_.end() ) { // Key already exists: Return it
			return *( i->second );
		} else { // Add a new key
			add( KeyC( id ) );
			return key( id );
		}
	}


	/// @brief Key with an identifier: Generate key if not present
	template< typename KeyC > // Concrete Key type to create
	inline
	static
	Key const &
	gen(
		std::string const & id,
		std::string const & identifier
	)
	{
		Map const & m_( m() );
		typename Map::const_iterator const i( m_.find( stripped_whitespace( id ) ) );
		if ( i != m_.end() ) { // Key already exists: Return it
			return *( i->second );
		} else { // Add a new key
			add( KeyC( id, identifier ) );
			return key( id );
		}
	}


	/// @brief Key with an identifier: Generate key if not present
	template< typename KeyC > // Concrete Key type to create
	inline
	static
	Key const &
	gen(
		std::string const & id,
		std::string const & identifier,
		std::string const & code
	)
	{
		Map const & m_( m() );
		typename Map::const_iterator const i( m_.find( stripped_whitespace( id ) ) );
		if ( i != m_.end() ) { // Key already exists: Return it
			return *( i->second );
		} else { // Add a new key
			add( KeyC( id, identifier, code ) );
			return key( id );
		}
	}


	/// @brief Remove a key
	inline
	static
	void
	remove( Key const & key )
	{
		// Remove all three identifiers from lookup map
		erase( key.id() );
		erase( key.identifier() );
		erase( key.code() );

		// Remove any others with same key from lookup map
		Map & m_( m() );
		for ( typename Map::iterator i = m_.begin(); i != m_.end(); ) {
			if ( i->second == key ) {
				m_.erase( i++ ); // Safe idiom for erasing while iterating through map
			} else {
				++i;
			}
		}

		// Remove key from collection
		keys().erase( key );
	}


	/// @brief Clear the collection and lookup map
	inline
	static
	void
	clear()
	{
		m().clear();
		keys().clear();
	}


public: // Properties


	/// @brief Size
	inline
	static
	Size
	size()
	{
		return keys().size();
	}


	/// @brief Empty?
	inline
	static
	bool
	empty()
	{
		return keys().empty();
	}


	/// @brief Has a key with an identifier?
	inline
	static
	bool
	has( std::string const & id )
	{
		Map const & m_( m() );
		return ( m_.find( stripped_whitespace( id ) ) != m_.end() );
	}


	/// @brief Has a key?
	inline
	static
	bool
	has( Key const & key )
	{
		Keys const & keys_( keys() );
		return ( keys_.find( key ) != keys_.end() );
	}


	/// @brief Key with an identifier
	inline
	static
	Key const &
	key( std::string const & id )
	{
		Map const & m_( m() );
		typename Map::const_iterator const i( m_.find( stripped_whitespace( id ) ) );
		debug_assert( i != m_.end() );
		return *( i->second );
	}


	/// @brief Number of keys with distinct indexes
	inline
	static
	typename Key::Size
	n_key()
	{
		return Key::n_key();
	}


public: // Iterators


	/// @brief Begin const iterator
	inline
	static
	ConstIterator
	begin()
	{
		return keys().begin();
	}


	/// @brief End const iterator
	inline
	static
	ConstIterator
	end()
	{
		return keys().end();
	}


private: // Methods


	/// @brief Insert a pair< identifier, KeyP > into the lookup map
	inline
	static
	void
	insert( std::string const & id, Key_ const & key )
	{
		if ( not_blank( id ) ) { // Don't insert if identifier is blank
			std::string const sid( stripped_whitespace( id ) );
			Map & m_( m() );
			debug_assert( ( m_.find( sid ) == m_.end() ) || ( *( m_.find( sid )->second ) == key ) ); // Catch repeat identifiers with unequal keys
			m_.insert( std::make_pair( sid, &key ) );
		}
	}


	/// @brief Insert a pair< identifier, KeyP > into the lookup map: Allow non-unique identifier
	inline
	static
	void
	insert_nonunique( std::string const & id, Key_ const & key )
	{
		if ( not_blank( id ) ) { // Don't insert if identifier is blank
			std::string const sid( stripped_whitespace( id ) );
			Map & m_( m() );
			m_.insert( std::make_pair( sid, &key ) );
		}
	}


	/// @brief Erase a pair< identifier, KeyP > from the lookup map
	inline
	static
	void
	erase( std::string const & id )
	{
		if ( not_blank( id ) ) { // No action if identifier is blank
			m().erase( stripped_whitespace( id ) );
		}
	}


	/// @brief Whitespace stripped tails copy of a string
	inline
	static
	std::string
	stripped_whitespace( std::string const & s )
	{
		using std::string;
		static string const WHITESPACE( " \t\0", 3 );
		if ( s.empty() ) return s; // Empty string
		string::size_type const b( s.find_first_not_of( WHITESPACE ) );
		string::size_type const e( s.find_last_not_of( WHITESPACE ) );
		if ( ( b == string::npos ) || ( e == string::npos ) ) { // Whitespace string
			return string(); // Return empty string
		} else {
			return s.substr( b, e - b + 1 ); // Trimmed
		}
	}


private: // Properties


	/// @brief Key collection
	inline
	static
	Keys &
	keys()
#ifdef __MINGW32__ // Work-around MinGW GCC 3.4.5 bug causing SIGSEGV with -O2 or -O3
	__attribute__ ((noinline))
#endif
	{
		static Keys keys_; // Function-local to support globals
		return keys_;
	}


	/// @brief Map from identifier to (non-owning) pointer to key
	inline
	static
	Map &
	m()
#ifdef __MINGW32__ // Work-around MinGW GCC 3.4.5 bug causing SIGSEGV with -O2 or -O3
	__attribute__ ((noinline))
#endif
	{
		static Map m_; // Function-local to support globals
		return m_;
	}


	/// @brief string is not blank?
	inline
	static
	bool
	not_blank( std::string const & s )
	{
		using std::string;
		static string const WHITESPACE( " \t\0", 3 );
		return ( s.find_first_not_of( WHITESPACE ) != string::npos );
	}


}; // KeyLookup


namespace lookup { // Lookup functors for key collection convenience


/// @brief Key lookup has query functor
template< typename K >
struct has
{


	/// @brief Default constructor
	inline
	has()
	= default;


	/// @brief Has a key with an identifier?
	inline
	bool
	operator ()( std::string const & id ) const
	{
		return KeyLookup< K >::has( id );
	}


	/// @brief Has a key?
	inline
	bool
	operator ()( K const & key ) const
	{
		return KeyLookup< K >::has( key );
	}


	/// @brief Has a key with the id of a key?
	inline
	bool
	operator ()( Key const & key_ ) const
	{
		return KeyLookup< K >::has( key_.id() );
	}


}; // has


/// @brief Key lookup key query functor
template< typename K >
struct key
{


	/// @brief Default constructor
	inline
	key()
	= default;


	/// @brief Key with an id
	inline
	K const &
	operator ()( std::string const & id ) const
	{
		return KeyLookup< K >::key( id );
	}


	/// @brief Key with the id of a key
	inline
	K const &
	operator ()( Key const & key_ ) const
	{
		return KeyLookup< K >::key( key_.id() );
	}


}; // key


/// @brief Key lookup/generator functor
template< typename K >
struct gen
{


	/// @brief Default constructor
	inline
	gen()
	= default;


	/// @brief Key with an id: Generate if not present
	template< typename KeyC > // Concrete Key type to create
	inline
	K const &
	key( std::string const & id ) const
	{
		return KeyLookup< K >::template gen< KeyC >( id );
	}


}; // gen


/// @brief Key lookup n_key query functor
template< typename K >
struct n_key
{


	/// @brief Default constructor
	inline
	n_key()
	= default;


	/// @brief Number of keys with distinct indexes
	inline
	typename K::Size
	operator ()() const
	{
		return K::n_key();
	}


}; // n_key


/// @brief Key lookup begin iterator functor
template< typename K >
struct begin
{


	// STL/boost style
	typedef  typename KeyLookup< K >::const_iterator  const_iterator;

	// Project style
	typedef  typename KeyLookup< K >::ConstIterator  ConstIterator;


	/// @brief Default constructor
	inline
	begin()
	= default;


	/// @brief Key collection begin iterator
	inline
	ConstIterator
	operator ()() const
	{
		return KeyLookup< K >::begin();
	}


}; // begin


/// @brief Key lookup end iterator functor
template< typename K >
struct end
{


	// STL/boost style
	typedef  typename KeyLookup< K >::const_iterator  const_iterator;

	// Project style
	typedef  typename KeyLookup< K >::ConstIterator  ConstIterator;


	/// @brief Default constructor
	inline
	end()
	= default;


	/// @brief Key collection end iterator
	inline
	ConstIterator
	operator ()() const
	{
		return KeyLookup< K >::end();
	}


}; // end


} // namespace lookup


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_KeyLookup_HH
