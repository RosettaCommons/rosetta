// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/StoredRotamerLibraryCreator.cc
/// @brief  Class for instantiating a SingleLigandRotamerLibrary from a StoredRotamerLibrarySpecification
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/rotamers/StoredRotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>

// Program header
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>
#include <core/chemical/ResidueType.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace rotamers {

static thread_local basic::Tracer TR("core.pack.rotamers.StoredRotamerLibraryCreator");

core::pack::rotamers::SingleResidueRotamerLibraryCOP
StoredRotamerLibraryCreator::create( core::chemical::ResidueType const & restype) const {
	using namespace core::chemical::rotamers;
	using namespace core::pack::dunbrack;

	RotamerLibrarySpecificationCOP libspec( restype.rotamer_library_specification() );
	// If the factory system is sound, these two checks should work.
	assert( libspec );
	StoredRotamerLibrarySpecificationCOP stored_libspec = utility::pointer::dynamic_pointer_cast< StoredRotamerLibrarySpecification const >(libspec);
	assert( stored_libspec );

	SingleLigandRotamerLibraryOP rotamer_lib( new SingleLigandRotamerLibrary() );
	rotamer_lib->init_from_vector( stored_libspec->coordinates() );
	return rotamer_lib;
}

std::string
StoredRotamerLibraryCreator::keyname() const {
	return core::chemical::rotamers::StoredRotamerLibrarySpecification::library_name();
}

} //namespace rotamers
} //namespace pack
} //namespace core

