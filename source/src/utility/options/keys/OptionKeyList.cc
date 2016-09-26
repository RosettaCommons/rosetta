// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/OptionKeyList.cc
/// @brief  Functions for adding option keys into an OptionKeyList
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <utility/options/keys/OptionKeyList.hh>

// Package headers
#include <utility/options/keys/all.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/keys/OptionKey.hh>

namespace utility {
namespace options {


/// @details For example if you have to add the option keys
/// packing::ex1, packing::ex2, packing::ex3 and packing::ex4
/// to a list of option keys named optlist, then you can write:
///
///     using namespace basic::options::OptionKeys::packing;
///     optlist + ex1 + ex2 + ex3 + ex4;
OptionKeyList &
operator + ( OptionKeyList & key_list, BooleanOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, IntegerOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, RealOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, StringOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, FileOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, PathOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, BooleanVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, IntegerVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, RealVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, StringVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, FileVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, PathVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, ResidueChainVectorOptionKey const & key )
{
	key_list.push_back( keys::VariantKey< OptionKey > ( key ) );
	return key_list;
}


} // namespace options
} // namespace utility

