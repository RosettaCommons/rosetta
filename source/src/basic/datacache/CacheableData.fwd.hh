// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableData.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_basic_datacache_CacheableData_fwd_hh
#define INCLUDED_basic_datacache_CacheableData_fwd_hh

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

class CacheableData;
typedef utility::pointer::shared_ptr< CacheableData > CacheableDataOP;
typedef utility::pointer::shared_ptr< CacheableData const > CacheableDataCOP;
typedef utility::pointer::weak_ptr< CacheableData > CacheableDataAP;
typedef utility::pointer::weak_ptr< CacheableData const > CacheableDataCAP;


} // namespace datacache
} // namespace basic


#endif /* INCLUDED_basic_datacache_CacheableData_FWD_HH */
