// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/WritableCacheableData.fwd.hh
/// @brief
/// @author Justin Porter


#ifndef INCLUDED_basic_datacache_WriteableCacheableData_fwd_hh
#define INCLUDED_basic_datacache_WriteableCacheableData_fwd_hh

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

class WriteableCacheableData;
typedef utility::pointer::shared_ptr< WriteableCacheableData > WriteableCacheableDataOP;
typedef utility::pointer::shared_ptr< WriteableCacheableData const > WriteeableCacheableDataCOP;
typedef utility::pointer::weak_ptr< WriteableCacheableData > WriteableCacheableDataAP;
typedef utility::pointer::weak_ptr< WriteableCacheableData const > WriteableCacheableDataCAP;


} // namespace datacache
} // namespace basic


#endif /* INCLUDED_basic_datacache_WritableCacheableData_FWD_HH */
