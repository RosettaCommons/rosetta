// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_basic_options_util_hh
#define INCLUDED_basic_options_util_hh

// C++ headers
#include <iosfwd>              // for string

// utility headers
#include <utility/vector1.hh>  // for vector1
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

namespace basic {
namespace options {

/// convenience functions

std::string
start_file();

utility::vector1< std::string >
start_files();

/// @brief Add a base-class OptionKey to an OptionCollection by trying to perform a dynamic
/// cast to each of the available option types until the correct derived class is found.
void
add_anonymous_option(
	utility::options::OptionCollection & options,
	utility::options::OptionKey const & key
);

/// @brief Create an OptionCollection that knows about only the subset of option keys
/// listed in the input OptionsKeys list. Load this new OptionCollection with the values
/// and the defaults that are stored in the global option collection.
utility::options::OptionCollectionOP
read_subset_of_global_option_collection(
	utility::options::OptionKeyList const & opt_keys
);

std::string
replace_option_namespace_colons_with_underscores( utility::options::OptionKey const & key );

} // namespace options
} // namespace basic


#endif // INCLUDED_basic_options_option_HH
