// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/CacheableAtomID_MapVector.cc
/// @brief  serialization routines for the CacheableAtomID_MapVector class
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/id/CacheableAtomID_MapVector.hh>

#ifdef SERIALIZATION
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

namespace core {
namespace id {

CacheableAtomID_MapVector::CacheableAtomID_MapVector() = default;

CacheableAtomID_MapVector::~CacheableAtomID_MapVector()= default;

basic::datacache::CacheableDataOP
CacheableAtomID_MapVector::clone() const {
	return basic::datacache::CacheableDataOP( new CacheableAtomID_MapVector(*this) );
}

}
}

#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::CacheableAtomID_MapVector::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( map_ ) ); // core::id::AtomID_Map<numeric::xyzVector<core::Real> >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::CacheableAtomID_MapVector::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( map_ ); // core::id::AtomID_Map<numeric::xyzVector<core::Real> >
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::CacheableAtomID_MapVector );
CEREAL_REGISTER_TYPE( core::id::CacheableAtomID_MapVector )

CEREAL_REGISTER_DYNAMIC_INIT( core_id_CacheableAtomID_MapVector )

#endif // SERIALIZATION
