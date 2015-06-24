// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/CenrotRotamerLibrarySpecification.cc
/// @brief  The CenrotRotamerLibrarySpecification class specifies building centroid rotamers.
/// @author Rocco Moretti (rmorettase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/CenrotRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/CenrotRotamerLibrarySpecificationCreator.hh>

// Utility headers
#include <utility/exit.hh>

namespace core {
namespace chemical {
namespace rotamers {

// Creator Functions

RotamerLibrarySpecificationOP
CenrotRotamerLibrarySpecificationCreator::create() const {
	return RotamerLibrarySpecificationOP( new CenrotRotamerLibrarySpecification );
}

RotamerLibrarySpecificationOP
CenrotRotamerLibrarySpecificationCreator::create( std::istream & input ) const {
	return RotamerLibrarySpecificationOP( new CenrotRotamerLibrarySpecification( input ) );
}

std::string
CenrotRotamerLibrarySpecificationCreator::keyname() const {
	return CenrotRotamerLibrarySpecification::library_name();
}

// Specification Functions

CenrotRotamerLibrarySpecification::CenrotRotamerLibrarySpecification():
	aa_( aa_unk )
{}

CenrotRotamerLibrarySpecification::CenrotRotamerLibrarySpecification( AA aa_setting )
{
	aa( aa_setting );
}

CenrotRotamerLibrarySpecification::CenrotRotamerLibrarySpecification( std::istream & input )
{
	std::string tag;
	input >> tag;
	if( ! input ) {
		utility_exit_with_message("Cannot find AA parameter on Cenrot rotamer input.");
	}
	aa( aa_from_name( tag ) );
}

void
CenrotRotamerLibrarySpecification::aa( AA aa_setting ) {
  if( aa_setting > num_canonical_aas ) {
    utility_exit_with_message("Cannot have a centroid rotamer library with a non-cannonical amino acid");
  }
	aa_ = aa_setting;
}
CenrotRotamerLibrarySpecification::~CenrotRotamerLibrarySpecification() {}

std::string
CenrotRotamerLibrarySpecification::keyname() const {
	return library_name();
}

std::string
CenrotRotamerLibrarySpecification::library_name() {
	return "CENROT";
}

} //namespace rotamers
} //namespace chemical
} //namespace core
