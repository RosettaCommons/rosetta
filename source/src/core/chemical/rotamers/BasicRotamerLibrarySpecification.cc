// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/BasicRotamerLibrarySpecification.cc
/// @brief  The BasicRotamerLibrarySpecification class says to build simple rotamer libraries
/// @author Rocco Moretti (rmorettase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/BasicRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/BasicRotamerLibrarySpecificationCreator.hh>

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

// Creator Functions

RotamerLibrarySpecificationOP
BasicRotamerLibrarySpecificationCreator::create() const {
	return RotamerLibrarySpecificationOP( new BasicRotamerLibrarySpecification );
}

RotamerLibrarySpecificationOP
BasicRotamerLibrarySpecificationCreator::create( std::istream & ) const {
	// We ignore any parameters
	return RotamerLibrarySpecificationOP( new BasicRotamerLibrarySpecification() );
}

std::string
BasicRotamerLibrarySpecificationCreator::keyname() const {
	return BasicRotamerLibrarySpecification::library_name();
}

// Specification Functions

BasicRotamerLibrarySpecification::BasicRotamerLibrarySpecification() {}

BasicRotamerLibrarySpecification::~BasicRotamerLibrarySpecification() {}

std::string
BasicRotamerLibrarySpecification::keyname() const {
	return library_name();
}

/// @details There's only a single BasicRotamerLibrary
std::string
BasicRotamerLibrarySpecification::cache_tag(ResidueType const &) const {
	return library_name();
}

std::string
BasicRotamerLibrarySpecification::library_name() {
	return "BASIC";
}

} //namespace rotamers
} //namespace chemical
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::rotamers::BasicRotamerLibrarySpecification::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::rotamers::RotamerLibrarySpecification >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::rotamers::BasicRotamerLibrarySpecification::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::rotamers::RotamerLibrarySpecification >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::rotamers::BasicRotamerLibrarySpecification );
CEREAL_REGISTER_TYPE( core::chemical::rotamers::BasicRotamerLibrarySpecification )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_rotamers_BasicRotamerLibrarySpecification )
#endif // SERIALIZATION
