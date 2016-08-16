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

StoredRotamerLibrarySpecification::StoredRotamerLibrarySpecification()
{}

StoredRotamerLibrarySpecification::StoredRotamerLibrarySpecification(std::istream & ) {
	utility_exit_with_message("Cannot currently instantiate a StoredRotamerLibrarySpecification from an input stream.");
}

StoredRotamerLibrarySpecification::~StoredRotamerLibrarySpecification() {}

/// @brief Add a particular rotamer to the list
void
StoredRotamerLibrarySpecification::add_rotamer( std::map< std::string, core::Vector > const & rotamer ) {
	coordinates_.push_back( rotamer );
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
