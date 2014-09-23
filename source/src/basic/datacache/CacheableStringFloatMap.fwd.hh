// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheableStringFloatMap.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_basic_datacache_CacheableStringFloatMap_fwd_hh
#define INCLUDED_basic_datacache_CacheableStringFloatMap_fwd_hh

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


class CacheableStringFloatMap;
typedef utility::pointer::shared_ptr< CacheableStringFloatMap > CacheableStringFloatMapOP;
typedef utility::pointer::shared_ptr< CacheableStringFloatMap const > CacheableStringFloatMapCOP;
typedef utility::pointer::weak_ptr< CacheableStringFloatMap > CacheableStringFloatMapAP;
typedef utility::pointer::weak_ptr< CacheableStringFloatMap const > CacheableStringFloatMapCAP;


} // namespace datacache
} // namespace basic

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#endif


#endif /* INCLUDED_basic_datacache_CacheableStringFloatMap_FWD_HH */
