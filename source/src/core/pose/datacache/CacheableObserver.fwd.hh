// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/datacache/CacheableObserver.fwd.hh
/// @brief  Forward declaration for CacheableObserver.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_datacache_CacheableObserver_fwd_hh
#define INCLUDED_core_pose_datacache_CacheableObserver_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace pose {
namespace datacache {


/// @brief forward declaration for CacheableObserver
class CacheableObserver;


/// @brief CacheableObserver const access pointer
typedef utility::pointer::shared_ptr< CacheableObserver const > CacheableObserverCAP;


/// @brief CacheableObserver access pointer
typedef utility::pointer::shared_ptr< CacheableObserver > CacheableObserverAP;


/// @brief CacheableObserver const owning pointer
typedef utility::pointer::shared_ptr< CacheableObserver const > CacheableObserverCOP;


/// @brief CacheableObserver owning pointer
typedef utility::pointer::shared_ptr< CacheableObserver > CacheableObserverOP;


} // namespace datacache
} // namespace pose
} // namespace core


#endif /* INCLUDED_core_pose_datacache_CacheableObserver_FWD_HH */
