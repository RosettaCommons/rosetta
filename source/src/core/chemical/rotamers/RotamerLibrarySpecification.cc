// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/RotamerLibrarySpecification.cc
/// @brief  The RotamerLibrarySpecification class tells how to build a rotamer library for a ResidueType
/// @author Rocco Moretti (rmorettase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Utility headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace rotamers {

RotamerLibrarySpecification::RotamerLibrarySpecification() = default;

RotamerLibrarySpecification::~RotamerLibrarySpecification() = default;

} //namespace rotamers
} //namespace chemical
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::rotamers::RotamerLibrarySpecification::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::rotamers::RotamerLibrarySpecification::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::rotamers::RotamerLibrarySpecification );
CEREAL_REGISTER_TYPE( core::chemical::rotamers::RotamerLibrarySpecification )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_rotamers_RotamerLibrarySpecification )
#endif // SERIALIZATION
