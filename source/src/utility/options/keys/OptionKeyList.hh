// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/OptionKeyList.hh
/// @brief  Functions for adding option keys into an OptionKeyList
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_options_keys_OptionKeyList_hh
#define INCLUDED_utility_options_keys_OptionKeyList_hh


// Unit headers
#include <utility/options/keys/OptionKeyList.fwd.hh>

// Package headers
#include <utility/options/keys/all.fwd.hh>

namespace utility {
namespace options {


/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, BooleanOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, IntegerOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, RealOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, StringOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, FileOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, PathOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, BooleanVectorOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, IntegerVectorOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, RealVectorOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, StringVectorOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, FileVectorOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, PathVectorOptionKey const & key );

/// @brief Convenience function for appending an option key to an option
/// key list (which itself is a list of pointers to OptionKeys, using
/// the VariantKey class) so that many option keys can be appended by
/// simply chaining a set of calls using "+".
OptionKeyList &
operator + ( OptionKeyList & key_list, ResidueChainVectorOptionKey const & key );


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_OptionKey_HH
