// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/RotamerLibrarySpecificationCreator.hh
/// @brief  Class for instantiating a particular RotamerLibrarySpecification
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_RotamerLibrarySpecificationCreator_HH
#define INCLUDED_core_chemical_rotamers_RotamerLibrarySpecificationCreator_HH

// Package headers
//#include <core/chemical/rotamers/RotamerLibrarySpecificationCreator.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.fwd.hh>

// Program header

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class RotamerLibrarySpecificationCreator : public utility::pointer::ReferenceCount {
public:
	virtual RotamerLibrarySpecificationOP
	create() const = 0;

	virtual RotamerLibrarySpecificationOP
	create( std::istream & ) const = 0;

	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< RotamerLibrarySpecificationCreator > RotamerLibrarySpecificationCreatorOP;
typedef utility::pointer::shared_ptr< RotamerLibrarySpecificationCreator const > RotamerLibrarySpecificationCreatorCOP;

} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
