// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/DunbrackRotamerLibrarySpecification.cc
/// @brief  The DunbrackRotamerLibrarySpecification class specifies building DunbrackRotamers.
/// @author Rocco Moretti (rmorettase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/DunbrackRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/DunbrackRotamerLibrarySpecificationCreator.hh>

// Utility headers
#include <utility/exit.hh>

namespace core {
namespace chemical {
namespace rotamers {

// Creator Functions

RotamerLibrarySpecificationOP
DunbrackRotamerLibrarySpecificationCreator::create() const {
	return RotamerLibrarySpecificationOP( new DunbrackRotamerLibrarySpecification );
}

RotamerLibrarySpecificationOP
DunbrackRotamerLibrarySpecificationCreator::create( std::istream & input ) const {
	return RotamerLibrarySpecificationOP( new DunbrackRotamerLibrarySpecification( input ) );
}

std::string
DunbrackRotamerLibrarySpecificationCreator::keyname() const {
	return DunbrackRotamerLibrarySpecification::library_name();
}

// Specification Functions

DunbrackRotamerLibrarySpecification::DunbrackRotamerLibrarySpecification():
		aa_(aa_unk)
{}

DunbrackRotamerLibrarySpecification::DunbrackRotamerLibrarySpecification( AA aa_setting )
{
	aa( aa_setting );
}

DunbrackRotamerLibrarySpecification::DunbrackRotamerLibrarySpecification( std::istream & input )
{
	std::string tag;
	input >> tag;
	if( ! input ) {
		utility_exit_with_message("Cannot find AA parameter on Dunbrack rotamer input line.");
	}
	aa( aa_from_name( tag ) );
}

DunbrackRotamerLibrarySpecification::~DunbrackRotamerLibrarySpecification() {}

void
DunbrackRotamerLibrarySpecification::aa( AA aa_setting ) {
  if( aa_setting > num_canonical_aas && ! core::chemical::is_canonical_D_aa(aa_setting) ) {
    utility_exit_with_message("Cannot have a Dunbrack rotamer library with a non-cannonical amino acid." );
  }
	aa_ = aa_setting;
}
std::string
DunbrackRotamerLibrarySpecification::keyname() const {
	return library_name();
}

std::string
DunbrackRotamerLibrarySpecification::library_name() {
	return "DUNBRACK";
}

} //namespace rotamers
} //namespace chemical
} //namespace core
