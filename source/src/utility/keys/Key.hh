// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key.hh
/// @brief  Hidden index key interface class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Interface provides common base class for covariant return in hierarchies
///  @li Derived classes specify the friend class(es) that can access the index
///  @li Index can be ignored if not used as key/index into a container
///  @li Can safely derive from concrete Key types as long as fields aren't added
///  @li Can derive privately from an Key to share the index set of another Key type
///      without allowing convertibility (for type safety)


#ifndef INCLUDED_utility_keys_Key_hh
#define INCLUDED_utility_keys_Key_hh


// Unit headers
#include <utility/keys/Key.fwd.hh>

// Package headers
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/KeyLess.fwd.hh>

// C++ headers
#include <cstddef>
#include <string>


namespace utility {
namespace keys {


/// @brief Hidden index key interface class
class Key
{


public: // Types


	// STL/boost style
	typedef  std::size_t  index_type;
	typedef  std::size_t  size_type;

	// Project style
	typedef  std::size_t  Index;
	typedef  std::size_t  Size;


private: // Friends


	template< typename O, typename S, typename C > friend class AutoKey;
	template< typename O, typename S, typename C > friend class UserKey;
	template< typename T, typename U > friend class KeyLess;
	template< typename T, typename U > friend class PointerKeyLess;


protected: // Creation


	/// @brief Default constructor
	inline
	Key()
	= default;


	/// @brief Copy constructor
	inline
	Key( Key const & )
	= default;


public: // Creation


	/// @brief Clone this
	virtual
	Key *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~Key()
	= default;


public: // Assignment


	/// @brief Copy assignment
	inline
	Key &
	operator =( Key const & key )
	{
		if ( this != &key ) {
			assign_Key( key );
		}
		return *this;
	}


protected: // Assignment


	/// @brief Key assignment
	virtual
	void
	assign_Key( Key const & key ) = 0;


public: // Properties


	/// @brief ID
	virtual
	std::string const &
	id() const = 0;


	/// @brief ID
	virtual
	std::string &
	id() = 0;


	/// @brief ID assignment
	virtual
	Key &
	id( std::string const & id_a ) = 0;


	/// @brief Identifier
	virtual
	std::string const &
	identifier() const = 0;


	/// @brief Identifier
	virtual
	std::string &
	identifier() = 0;


	/// @brief Identifier assignment
	virtual
	Key &
	identifier( std::string const & identifier_a ) = 0;


	/// @brief Code
	virtual
	std::string const &
	code() const = 0;


	/// @brief Code
	virtual
	std::string &
	code() = 0;


	/// @brief Code assignment
	virtual
	Key &
	code( std::string const & code_a ) = 0;


	/// @brief Index
	/// @note  Only for use as an optimization: DO NOT WRITE CODE DEPENDING ON THE SPECIFIC INDEX VALUE!
	virtual
	Index
	private_index() const = 0;


public: // Comparison


	/// @brief Key == Key
	/// @note  TYpe and index-based polymorphic equality
	friend
	inline
	bool
	operator ==( Key const & a, Key const & b )
	{
		return a.equals( b );
	}


	/// @brief Key != Key
	/// @note  TYpe and index-based polymorphic equality
	friend
	inline
	bool
	operator !=( Key const & a, Key const & b )
	{
		return ! a.equals( b );
	}


	/// @brief Key < Key
	/// @note  Index-based ordering
	/// @note  Needed for use as key in associative containers
	friend
	inline
	bool
	operator <( Key const & a, Key const & b )
	{
		return a.less_than( b );
	}


	/// @brief Key <= Key
	/// @note  Index-based ordering
	/// @note  Needed for use as key in associative containers
	friend
	inline
	bool
	operator <=( Key const & a, Key const & b )
	{
		return ( ( a.less_than( b ) ) || ( a.equals( b ) ) );
	}


	/// @brief Key >= Key
	/// @note  Index-based ordering
	/// @note  Needed for use as key in associative containers
	friend
	inline
	bool
	operator >=( Key const & a, Key const & b )
	{
		return ( ( b.less_than( a ) ) || ( a.equals( b ) ) );
	}


	/// @brief Key > Key
	/// @note  Index-based ordering
	/// @note  Needed for use as key in associative containers
	friend
	inline
	bool
	operator >( Key const & a, Key const & b )
	{
		return b.less_than( a );
	}


	/// @brief Are Keys of Comparable Types?
	friend
	inline
	bool
	comparable( Key const & a, Key const & b )
	{
		return a.comparable( b );
	}


protected: // Conversion


	/// @brief Index value conversion
	/// @note  A pure virtual version of this slows down key lookup operations
	///        because it prevents inlining for derived key types
	inline
	operator Index() const
	{
		return index();
	}


protected: // Properties


	/// @brief Index
	virtual
	Index
	index() const = 0;


	/// @brief Index
	virtual
	Index &
	index() = 0;


	/// @brief Index assignment
	virtual
	Key &
	index( Index const index_a ) = 0;


protected: // Comparison


	/// @brief Equal to a Key?
	virtual
	bool
	equals( Key const & key ) const = 0;


	/// @brief Less than a Key?
	virtual
	bool
	less_than( Key const & key ) const = 0;


	/// @brief Comparable to a Key?
	virtual
	bool
	comparable( Key const & key ) const = 0;


}; // Key


// Friend function namespace declarations


/// @brief Key == Key
bool
operator ==( Key const & a, Key const & b );


/// @brief Key != Key
bool
operator !=( Key const & a, Key const & b );


/// @brief Key < Key
bool
operator <( Key const & a, Key const & b );


/// @brief Key <= Key
bool
operator <=( Key const & a, Key const & b );


/// @brief Key >= Key
bool
operator >=( Key const & a, Key const & b );


/// @brief Key > Key
bool
operator >( Key const & a, Key const & b );


/// @brief Are Keys of Comparable Types?
bool
comparable( Key const & a, Key const & b );


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_Key_HH
