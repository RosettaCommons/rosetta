// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/BasicDataCache.fwd.hh
/// @brief  forward declaration for BasicDataCache
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_basic_datacache_BasicDataCache_fwd_hh
#define INCLUDED_basic_datacache_BasicDataCache_fwd_hh

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


/// @brief forward declaration for BasicDataCache
class BasicDataCache;


/// @brief BasicDataCache owning pointer
typedef utility::pointer::shared_ptr< BasicDataCache > BasicDataCacheOP;


/// @brief BasicDataCache const owning pointer
typedef utility::pointer::shared_ptr< BasicDataCache const > BasicDataCacheCOP;


/// @brief BasicDataCache access pointer
typedef utility::pointer::weak_ptr< BasicDataCache > BasicDataCacheAP;


/// @brief BasicDataCache const access pointer
typedef utility::pointer::weak_ptr< BasicDataCache const > BasicDataCacheCAP;


} // namespace datacache
} // namespace basic


#endif /* INCLUDED_basic_datacache_BasicDataCache_FWD_HH */
