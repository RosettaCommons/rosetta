// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/OptionKey.hh
/// @brief  Abstract automatic hidden index key for options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_OptionKey_hh
#define INCLUDED_utility_options_keys_OptionKey_hh


// Unit headers
#include <utility/options/keys/OptionKey.fwd.hh>

// Package headers
#include <utility/options/Option.fwd.hh>

// Project headers
#include <utility/keys/AutoKey.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/VariantKey.hh>


namespace utility {
namespace options {


/// @brief Abstract automatic hidden index key for options
class OptionKey :
	public utility::keys::AutoKey< Option >
{


private: // Types


	typedef  utility::keys::AutoKey< Option >  Super;


protected: // Types


	typedef  utility::keys::Key  Key;


public: // Types


	typedef  utility::keys::KeyLookup< OptionKey >  Lookup;


private: // Friends


#if !(defined _MSC_VER) || (defined __INTEL_COMPILER) // Visual C++ 2005 bug work-around
	template< typename K, typename T > friend class utility::keys::SmallKeyVector;
#endif


protected: // Creation


	/// @brief Default constructor
	inline
	OptionKey()
	{}


	/// @brief Copy constructor
	inline
	OptionKey( OptionKey const & key ) :
		Super( key )
	{}


	/// @brief Copy + identifier constructor
	inline
	OptionKey(
		OptionKey const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{}


	/// @brief Key constructor
	inline
	explicit
	OptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	OptionKey(
		Key const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{}


	/// @brief Identifier constructor
	inline
	explicit
	OptionKey(
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{}


public: // Creation


	/// @brief Clone this
	virtual
	OptionKey *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~OptionKey() {}


public: // Assignment


	/// @brief Key assignment
	inline
	OptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


public: // Properties


	/// @brief Scalar option key?
	virtual
	bool
	scalar() const = 0;


	/// @brief Vector option key?
	virtual
	bool
	vector() const = 0;


}; // OptionKey

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
///
/// @details For example if you have to add the option keys
/// packing::ex1, packing::ex2, packing::ex3 and packing::ex4
/// to a list of option keys named optlist, then you can write:
///
///     using namespace basic::options::OptionKeys::packing;
///     optlist + ex1 + ex2 + ex3 + ex4;
inline
OptionKeyList &
operator + ( OptionKeyList & key_list, OptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_OptionKey_HH
