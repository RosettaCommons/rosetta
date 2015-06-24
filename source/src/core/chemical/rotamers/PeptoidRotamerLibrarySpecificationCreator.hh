// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/PeptoidRotamerLibrarySpecificationCreator.hh
/// @brief  Class for instantiating a particular PeptoidRotamerLibrarySpecification
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_PeptoidRotamerLibrarySpecificationCreator_HH
#define INCLUDED_core_chemical_rotamers_PeptoidRotamerLibrarySpecificationCreator_HH

// Package headers
#include <core/chemical/rotamers/RotamerLibrarySpecificationCreator.hh>
#include <core/chemical/rotamers/PeptoidRotamerLibrarySpecification.fwd.hh>


// Program header

// Utility headers

// C++ headers
#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class PeptoidRotamerLibrarySpecificationCreator : public RotamerLibrarySpecificationCreator {
public:
	virtual RotamerLibrarySpecificationOP
	create() const;

	virtual RotamerLibrarySpecificationOP
	create( std::istream & ) const;

	virtual std::string keyname() const;

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
