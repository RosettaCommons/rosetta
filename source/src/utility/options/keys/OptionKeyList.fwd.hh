// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/OptionKey.fwd.hh
/// @brief  utility::options::OptionKeyList forward declarations
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_options_keys_OptionKeyList_fwd_hh
#define INCLUDED_utility_options_keys_OptionKeyList_fwd_hh

#include <list>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>

namespace utility {
namespace options {

typedef std::list< keys::VariantKey< OptionKey > > OptionKeyList;

} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_OptionKeyList_FWD_HH
