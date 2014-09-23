// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheableAtomID_MapVector.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_id_CacheableAtomID_MapVector_fwd_hh
#define INCLUDED_core_id_CacheableAtomID_MapVector_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace id {


class CacheableAtomID_MapVector;
typedef utility::pointer::shared_ptr< CacheableAtomID_MapVector > CacheableAtomID_MapVectorOP;
typedef utility::pointer::shared_ptr< CacheableAtomID_MapVector const > CacheableAtomID_MapVectorCOP;
typedef utility::pointer::weak_ptr< CacheableAtomID_MapVector > CacheableAtomID_MapVectorAP;
typedef utility::pointer::weak_ptr< CacheableAtomID_MapVector const > CacheableAtomID_MapVectorCAP;


} // namespace id
} // namespace core


#endif /* INCLUDED_core_util_datacache_CacheableAtomID_MapVector_FWD_HH */
