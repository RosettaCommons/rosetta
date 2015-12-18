// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/datacache/CacheablePoseRawPtr.cc
/// @brief
/// @author Phil Bradley

// Unit headers
#include <core/pose/datacache/CacheablePoseRawPtr.hh>

// Package headers
#include <core/pose/Pose.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {

CacheablePoseRawPtr::CacheablePoseRawPtr( core::pose::PoseOP pose )
: CacheableData(), pose_(pose)
{}

CacheablePoseRawPtr::~CacheablePoseRawPtr(){}

basic::datacache::CacheableDataOP
CacheablePoseRawPtr::clone() const { return basic::datacache::CacheableDataOP( new CacheablePoseRawPtr(*this) ); }

core::pose::PoseOP CacheablePoseRawPtr::pose() { return pose_; }

}
}
}


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::pose::datacache::CacheablePoseRawPtr::CacheablePoseRawPtr() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::datacache::CacheablePoseRawPtr::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( pose_ ) ); // core::pose::PoseOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::datacache::CacheablePoseRawPtr::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( pose_ ); // core::pose::PoseOP
}
SAVE_AND_LOAD_SERIALIZABLE( core::pose::datacache::CacheablePoseRawPtr );
CEREAL_REGISTER_TYPE( core::pose::datacache::CacheablePoseRawPtr )


CEREAL_REGISTER_DYNAMIC_INIT( core_pose_datacache_CacheablePoseRawPtr )
#endif // SERIALIZATION



