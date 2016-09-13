// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/LongRangeEnergyContainer.cc
/// @brief  A container interface for storing and scoring long range energies
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/LREnergyContainer.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

ResidueNeighborIterator::~ResidueNeighborIterator() = default;

ResidueNeighborConstIterator::~ResidueNeighborConstIterator() = default;

LREnergyContainer::~LREnergyContainer() = default;

}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::LREnergyContainer::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::LREnergyContainer::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::LREnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::LREnergyContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_LREnergyContainer )
#endif // SERIALIZATION
