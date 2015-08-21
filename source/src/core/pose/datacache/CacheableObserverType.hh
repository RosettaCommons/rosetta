// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/datacache/CacheableObserverType.hh
/// @brief  enum indexing the observer types stored in a Pose's internal ObserverCache
/// @author

#ifndef INCLUDED_core_pose_datacache_CacheableObserverType_hh
#define INCLUDED_core_pose_datacache_CacheableObserverType_hh


namespace core {
namespace pose {
namespace datacache {


// hold the enum within a descriptive namespace to avoid name collisions
class CacheableObserverType
{
public:

	/// @brief enum indexing the data types stored in a Pose's internal DataCache
	enum Enum {
		// Remember to set the first enum in the list to 1!
		LENGTH_EVENT_COLLECTOR = 1,
		SPECIAL_SEGMENTS_OBSERVER,
		ENZDES_OBSERVER,
		// *** IMPORTANT ***
		// The 'num_cacheable_data_types' below must be the last enum, and must
		// always be set equal to the (last-1) enum.  If you append a new enum
		// to the list, remember to change the value below!
		num_cacheable_data_types = ENZDES_OBSERVER
	};

}; // class CacheableObserverType


} // namespace datacache
} // namespace pose
} // namespace core


#endif /* INCLUDED_core_pose_datacache_CacheableObserverType_HH */
