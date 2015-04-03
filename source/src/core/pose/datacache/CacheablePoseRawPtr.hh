// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheablePoseRawPtr.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_core_pose_datacache_CacheablePoseRawPtr_hh
#define INCLUDED_core_pose_datacache_CacheablePoseRawPtr_hh

// unit headers
#include <core/pose/datacache/CacheablePoseRawPtr.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>

// project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace datacache {


/// @brief pose *
class CacheablePoseRawPtr : public basic::datacache::CacheableData
{
public:
	CacheablePoseRawPtr( core::pose::Pose * pose ) : CacheableData(), pose_(pose) {}
	virtual ~CacheablePoseRawPtr(){};
	virtual basic::datacache::CacheableDataOP clone() const { return basic::datacache::CacheableDataOP( new CacheablePoseRawPtr(*this) ); }
	virtual core::pose::Pose * pose() const { return pose_; }
private:
	core::pose::Pose * pose_;
};


} // namespace datacache
} // namespace util
} // namespace core

#endif /* INCLUDED_core_pose_datacache_CacheablePoseRawPtr_HH */
