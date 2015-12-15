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


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {


/// @brief Holds a smart pointer (no longer a raw pointer) to a Pose so that it
/// can be tucked inside another Pose.
class CacheablePoseRawPtr : public basic::datacache::CacheableData
{
public:
	CacheablePoseRawPtr( core::pose::PoseOP pose );
	virtual ~CacheablePoseRawPtr();
	virtual basic::datacache::CacheableDataOP clone() const;
	virtual core::pose::PoseOP pose();

private:
	core::pose::PoseOP pose_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CacheablePoseRawPtr();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace util
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_datacache_CacheablePoseRawPtr )
#endif // SERIALIZATION


#endif /* INCLUDED_core_pose_datacache_CacheablePoseRawPtr_HH */
