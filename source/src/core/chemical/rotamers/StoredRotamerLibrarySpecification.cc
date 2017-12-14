// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/StoredRotamerLibrarySpecification.cc
/// @brief  The StoredRotamerLibrarySpecification class specifies building PDBRotamers.
/// @author Rocco Moretti (rmorettase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecificationCreator.hh>

#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/xyzVector.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace rotamers {

// Creator Functions

RotamerLibrarySpecificationOP
StoredRotamerLibrarySpecificationCreator::create() const {
	return RotamerLibrarySpecificationOP( new StoredRotamerLibrarySpecification );
}

RotamerLibrarySpecificationOP
StoredRotamerLibrarySpecificationCreator::create( std::istream & input ) const {
	return RotamerLibrarySpecificationOP( new StoredRotamerLibrarySpecification( input ) );
}

std::string
StoredRotamerLibrarySpecificationCreator::keyname() const {
	return StoredRotamerLibrarySpecification::library_name();
}

// Specification Functions

StoredRotamerLibrarySpecification::StoredRotamerLibrarySpecification() = default;

//NOTE You can use SingleLigandRotamerLibrary's utility function rotamer_information_from_PDB_stream to get objects you can pass to add_rotamers and set_reference_energy; see this class's unit test for an example.  The utility function can't be called here due to library level rules.
StoredRotamerLibrarySpecification::StoredRotamerLibrarySpecification(std::istream & ) {
	utility_exit_with_message("Cannot currently instantiate a StoredRotamerLibrarySpecification from an input stream.");
}

StoredRotamerLibrarySpecification::~StoredRotamerLibrarySpecification() = default;

/// @brief Add a particular rotamer to the list
void
StoredRotamerLibrarySpecification::add_rotamer( std::map< std::string, core::Vector > const & rotamer ) {
	coordinates_.push_back( rotamer );
}

/// @brief Add a vector of rotamers to the list
void
StoredRotamerLibrarySpecification::add_rotamers( utility::vector1< std::map< std::string, core::Vector > > const & rotamers ) {
	for ( auto const & rotamer : rotamers ) coordinates_.push_back( rotamer );
}

std::string
StoredRotamerLibrarySpecification::keyname() const {
	return library_name();
}

std::string
StoredRotamerLibrarySpecification::library_name() {
	return "STORED";
}

} //namespace rotamers
} //namespace chemical
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::rotamers::StoredRotamerLibrarySpecification::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::rotamers::RotamerLibrarySpecification >( this ) );
	arc( CEREAL_NVP( coordinates_ ) ); // utility::vector1<std::map<std::string, core::Vector> >
	arc( CEREAL_NVP( ref_energy_  ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::rotamers::StoredRotamerLibrarySpecification::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::rotamers::RotamerLibrarySpecification >( this ) );
	arc( coordinates_ ); // utility::vector1<std::map<std::string, core::Vector> >
	arc( ref_energy_  ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::rotamers::StoredRotamerLibrarySpecification );
CEREAL_REGISTER_TYPE( core::chemical::rotamers::StoredRotamerLibrarySpecification )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_rotamers_StoredRotamerLibrarySpecification )
#endif // SERIALIZATION
