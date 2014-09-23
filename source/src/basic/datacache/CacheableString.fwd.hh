// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheableString.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_basic_datacache_CacheableString_fwd_hh
#define INCLUDED_basic_datacache_CacheableString_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/down_cast.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>


namespace basic {
namespace datacache {


class CacheableString;
typedef utility::pointer::shared_ptr< CacheableString > CacheableStringOP;
typedef utility::pointer::shared_ptr< CacheableString const > CacheableStringCOP;
typedef utility::pointer::weak_ptr< CacheableString > CacheableStringAP;
typedef utility::pointer::weak_ptr< CacheableString const > CacheableStringCAP;


} // namespace datacache
} // namespace basic


#endif /* INCLUDED_basic_datacache_CacheableString_FWD_HH */
