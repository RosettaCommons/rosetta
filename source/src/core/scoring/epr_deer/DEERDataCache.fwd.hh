// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/epr_deer/DEERDataCache.cc
/// @brief  Contains data for DEER energy method, stored as DEERData objects
/// @details To prevent the energy method from reading from the command line every scoring
///      round and parsing the input file, this method stores a list of DEER decay data.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_DEERDataCache_fwd_hh
#define INCLUDED_core_scoring_epr_deer_DEERDataCache_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace epr_deer {

class DEERDataCache;

typedef utility::pointer::shared_ptr< DEERDataCache > DEERDataCacheOP;
typedef utility::pointer::shared_ptr< DEERDataCache const > DEERDataCacheCOP;

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
