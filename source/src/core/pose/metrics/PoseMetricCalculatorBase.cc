// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/PoseMetricCalculatorBase.cc
/// @brief  Auto-generated serialization template functions
/// @author Andrew Leaver-Fay (aleavefay@gmail.com)

// Unit headers
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::EnergyDependentCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::PoseMetricCalculator >( this ) );
	arc( CEREAL_NVP( energies_are_outdated_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::EnergyDependentCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::PoseMetricCalculator >( this ) );
	arc( energies_are_outdated_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::EnergyDependentCalculator );
CEREAL_REGISTER_TYPE( core::pose::metrics::EnergyDependentCalculator )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::StructureDependentCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::PoseMetricCalculator >( this ) );
	arc( CEREAL_NVP( structure_is_outdated_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::StructureDependentCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::PoseMetricCalculator >( this ) );
	arc( structure_is_outdated_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::StructureDependentCalculator );
CEREAL_REGISTER_TYPE( core::pose::metrics::StructureDependentCalculator )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::PoseMetricCalculator::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::PoseMetricCalculator::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::PoseMetricCalculator );
CEREAL_REGISTER_TYPE( core::pose::metrics::PoseMetricCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_metrics_PoseMetricCalculatorBase )
#endif // SERIALIZATION
