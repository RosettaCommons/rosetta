// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/CstMinimizationData.cc
/// @brief  A cacheable data wrapper for ConstraintsOPs for use in minimization
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/constraints/CstMinimizationData.hh>

// Package headers
#include <core/scoring/constraints/Constraints.hh>
#include <utility>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace constraints {

CstMinimizationData::CstMinimizationData() = default;

CstMinimizationData::CstMinimizationData( ConstraintsOP constraints ) : constraints_(std::move( constraints )) {}

CstMinimizationData::~CstMinimizationData() = default;

CstMinimizationData::CacheableDataOP
CstMinimizationData::clone() const { return CacheableDataOP( new CstMinimizationData( *this ) ); }

void CstMinimizationData::set_constraints( ConstraintsOP constraints ) { constraints_ = constraints; }

}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::CstMinimizationData::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( constraints_ ) ); // ConstraintsOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::CstMinimizationData::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( constraints_ ); // ConstraintsOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::CstMinimizationData );
CEREAL_REGISTER_TYPE( core::scoring::constraints::CstMinimizationData )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_CstMinimizationData )
#endif // SERIALIZATION
