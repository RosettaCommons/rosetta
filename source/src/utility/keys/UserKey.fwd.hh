// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/UserKey.fwd.hh
/// @brief  utility::keys::UserKey forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_keys_UserKey_fwd_hh
#define INCLUDED_utility_keys_UserKey_fwd_hh


// Package headers
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/NoClient.fwd.hh>


namespace utility {
namespace keys {


// Forward
template< typename O, typename S = Key, typename C = NoClient > class UserKey;


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_UserKey_FWD_HH
