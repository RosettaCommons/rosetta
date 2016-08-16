// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/datacache/cacheable_observers.fwd.hh
/// @brief  Forward declarations for a bunch of CacheableObserver implementations.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_datacache_cacheable_observers_fwd_hh
#define INCLUDED_core_pose_datacache_cacheable_observers_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace pose {
namespace datacache {


/// @brief forward declaration for LengthEventCollector
class LengthEventCollector;
typedef utility::pointer::shared_ptr< LengthEventCollector const > LengthEventCollectorCAP;
typedef utility::pointer::shared_ptr< LengthEventCollector > LengthEventCollectorAP;
typedef utility::pointer::shared_ptr< LengthEventCollector const > LengthEventCollectorCOP;
typedef utility::pointer::shared_ptr< LengthEventCollector > LengthEventCollectorOP;

class SpecialSegmentsObserver;
typedef utility::pointer::shared_ptr< SpecialSegmentsObserver const > SpecialSegmentsObserverCOP;
typedef utility::pointer::shared_ptr< SpecialSegmentsObserver > SpecialSegmentsObserverOP;

} // namespace datacache
} // namespace pose
} // namespace core


#endif /* INCLUDED_core_util_datacache_cacheable_observers_FWD_HH */
