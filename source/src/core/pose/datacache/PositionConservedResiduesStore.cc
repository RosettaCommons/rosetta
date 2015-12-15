// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/datacache/PositionConservedResiduesStore.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <core/pose/datacache/PositionConservedResiduesStore.hh>

// Project headers
#include <core/types.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/unordered_map.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace datacache {

PositionConservedResiduesStore::PositionConservedResiduesStore() {}

basic::datacache::CacheableDataOP PositionConservedResiduesStore::clone() const {
	return basic::datacache::CacheableDataOP( new PositionConservedResiduesStore(*this) );
}

bool PositionConservedResiduesStore::is_conserved(core::Size residue) const {
	return (conservation_.find(residue) != conservation_.end()) ? conservation_[residue] : false;
}

void PositionConservedResiduesStore::set_conserved(core::Size residue, bool conserved) {
	conservation_[residue] = conserved;
}

}  // namespace datacache
}  // namespace pose
}  // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::datacache::PositionConservedResiduesStore::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( conservation_ ) ); // ConservationMap
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::datacache::PositionConservedResiduesStore::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( conservation_ ); // ConservationMap
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::datacache::PositionConservedResiduesStore );
CEREAL_REGISTER_TYPE( core::pose::datacache::PositionConservedResiduesStore )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_datacache_PositionConservedResiduesStore )
#endif // SERIALIZATION
