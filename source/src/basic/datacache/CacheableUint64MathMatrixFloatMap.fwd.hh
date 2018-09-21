// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableUint64MathMatrixFloatMap.fwd.hh
/// @brief
/// @author Brian Coventry


#ifndef INCLUDED_basic_datacache_CacheableUint64MathMatrixFloatMap_fwd_hh
#define INCLUDED_basic_datacache_CacheableUint64MathMatrixFloatMap_fwd_hh

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


class CacheableUint64MathMatrixFloatMap;
typedef utility::pointer::shared_ptr< CacheableUint64MathMatrixFloatMap > CacheableUint64MathMatrixFloatMapOP;
typedef utility::pointer::shared_ptr< CacheableUint64MathMatrixFloatMap const > CacheableUint64MathMatrixFloatMapCOP;
typedef utility::pointer::weak_ptr< CacheableUint64MathMatrixFloatMap > CacheableUint64MathMatrixFloatMapAP;
typedef utility::pointer::weak_ptr< CacheableUint64MathMatrixFloatMap const > CacheableUint64MathMatrixFloatMapCAP;


} // namespace datacache
} // namespace basic


#endif /* INCLUDED_basic_datacache_CacheableUint64MathMatrixFloatMap_FWD_HH */
