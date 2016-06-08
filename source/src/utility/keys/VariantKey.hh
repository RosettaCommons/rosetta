// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/keys/VariantKey.hh
/// @brief  Variant key class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_keys_VariantKey_hh
#define INCLUDED_utility_keys_VariantKey_hh


// Unit headers
#include <utility/keys/VariantKey.fwd.hh>

// C++ headers
#include <utility/assert.hh>
#include <string>


namespace utility {
namespace keys {


/// @brief Variant key class
template< typename K >
class VariantKey
{


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  typename K::index_type  index_type;

	// Project style
	typedef  K  Key;
	typedef  typename K::Index  Index;


public: // Creation


	/// @brief Default constructor
	inline
	VariantKey() :
		key_p_( 0 )
	{}


	/// @brief Copy constructor
	inline
	VariantKey( VariantKey const & var ) :
		key_p_( var.key_p_ ? var.key_p_->clone() : 0 )
	{}


	/// @brief Key constructor
	inline
	VariantKey( Key const & key_a ) :
		key_p_( key_a.clone() )
	{}


	/// @brief Destructor
	inline
	~VariantKey()
	{
		delete key_p_;
	}


public: // Assignment


	/// @brief Copy assignment
	inline
	VariantKey &
	operator =( VariantKey const & var )
	{
		if ( this != &var ) {
			delete key_p_; key_p_ = ( var.key_p_ ? var.key_p_->clone() : 0 );
		}
		return *this;
	}


public: // Conversion


	/// @brief Key conversion
	inline
	operator Key const &() const
	{
		debug_assert( key_p_ );
		return *key_p_;
	}


	/// @brief Key conversion
	inline
	operator Key &()
	{
		debug_assert( key_p_ );
		return *key_p_;
	}


public: // Properties


	/// @brief ID
	inline
	std::string const &
	id() const
	{
		debug_assert( key_p_ );
		return key_p_->id();
	}


	/// @brief ID
	inline
	std::string &
	id()
	{
		debug_assert( key_p_ );
		return key_p_->id();
	}


	/// @brief ID assignment
	inline
	VariantKey &
	id( std::string const & id_a )
	{
		debug_assert( key_p_ );
		key_p_->id( id_a );
		return *this;
	}


	/// @brief Identifier
	inline
	std::string const &
	identifier() const
	{
		debug_assert( key_p_ );
		return key_p_->identifier();
	}


	/// @brief Identifier
	inline
	std::string &
	identifier()
	{
		debug_assert( key_p_ );
		return key_p_->identifier();
	}


	/// @brief Identifier assignment
	inline
	VariantKey &
	identifier( std::string const & identifier_a )
	{
		debug_assert( key_p_ );
		key_p_->identifier( identifier_a );
		return *this;
	}


	/// @brief Code
	inline
	std::string const &
	code() const
	{
		debug_assert( key_p_ );
		return key_p_->code();
	}


	/// @brief Code
	inline
	std::string &
	code()
	{
		debug_assert( key_p_ );
		return key_p_->code();
	}


	/// @brief Code assignment
	inline
	VariantKey &
	code( std::string const & code_a )
	{
		debug_assert( key_p_ );
		key_p_->code( code_a );
		return *this;
	}


	/// @brief Index
	/// @note  Only for use as an optimization: DO NOT WRITE CODE DEPENDING ON THE SPECIFIC INDEX VALUE!
	inline
	Index
	private_index() const
	{
		debug_assert( key_p_ );
		return key_p_->private_index();
	}


	/// @brief Key
	inline
	Key const &
	operator ()() const
	{
		debug_assert( key_p_ );
		return *key_p_;
	}


	/// @brief Key
	inline
	Key &
	operator ()()
	{
		debug_assert( key_p_ );
		return *key_p_;
	}


public: // Comparison


	/// @brief VariantKey == VariantKey
	friend
	inline
	bool
	operator ==( VariantKey const & a, VariantKey const & b )
	{
		return ( *a.key_p_ == *b.key_p_ );
	}


	/// @brief VariantKey != VariantKey
	friend
	inline
	bool
	operator !=( VariantKey const & a, VariantKey const & b )
	{
		return ( *a.key_p_ != *b.key_p_ );
	}


	/// @brief VariantKey < VariantKey
	friend
	inline
	bool
	operator <( VariantKey const & a, VariantKey const & b )
	{
		return ( *a.key_p_ < *b.key_p_ );
	}


private: // Conversion

	/* // Unused
	/// @brief Index value conversion
	/// @note  A pure virtual version of this slows down key lookup operations
	///        because it prevents inlining for derived key types
	inline
	operator Index() const
	{
		debug_assert( key_p_ );
		return key_p_->index();
	}
	*/

private: // Fields


	/// @brief Pointer to key
	Key * key_p_;


}; // VariantKey


// Friend function namespace declarations


/// @brief VariantKey == VariantKey
template< typename K >
bool
operator ==( VariantKey< K > const & a, VariantKey< K > const & b );


/// @brief VariantKey != VariantKey
template< typename K >
bool
operator !=( VariantKey< K > const & a, VariantKey< K > const & b );


/// @brief VariantKey < VariantKey
template< typename K >
bool
operator <( VariantKey< K > const & a, VariantKey< K > const & b );


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_VariantKey_HH
