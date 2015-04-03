// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/BasicDataCache.cc
/// @brief  A DataCache storing objects derived from
///         basic::datacache::CacheableData.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <basic/datacache/BasicDataCache.hh>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <vector>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>


namespace basic {
namespace datacache {


/// @brief default constructor
BasicDataCache::BasicDataCache() :
	Super()
{}


/// @brief size constructor
/// @param[in] n_types The number of slots for this DataCache.
BasicDataCache::BasicDataCache( size_t const n_slots ) :
	Super( n_slots )
{}


/// @brief copy constructor
BasicDataCache::BasicDataCache( BasicDataCache const & rval ) :
	Super( rval )
{}


/// @brief default destructor
BasicDataCache::~BasicDataCache() {}


/// @brief copy assignment
BasicDataCache & BasicDataCache::operator =( BasicDataCache const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
	}
	return *this;
}


} // namespace datacache
} // namespace basic
