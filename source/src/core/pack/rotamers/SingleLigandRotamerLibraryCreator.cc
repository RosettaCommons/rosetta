// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/SingleLigandRotamerLibraryCreator.cc
/// @brief  Class for instantiating a particular SingleResidueRotamerLibrary
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/rotamers/SingleLigandRotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>

// Program header
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/PDBRotamerLibrarySpecification.hh>
#include <core/chemical/ResidueType.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace rotamers {

static thread_local basic::Tracer TR("core.pack.rotamers.SingleLigandRotamerLibraryCreator");

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleLigandRotamerLibraryCreator::create( core::chemical::ResidueType const & restype) const {
	using namespace core::chemical::rotamers;
	using namespace core::pack::dunbrack;

	RotamerLibrarySpecificationCOP libspec( restype.rotamer_library_specification() );
	// If the factory system is sound, these two checks should work.
	assert( libspec );
	PDBRotamerLibrarySpecificationCOP pdb_libspec = utility::pointer::dynamic_pointer_cast< PDBRotamerLibrarySpecification const >(libspec);
	assert( pdb_libspec );

	std::string const & rotamers_file(pdb_libspec->pdb_rotamers_file());

	if ( rotamers_file.size() == 0 ) {
		utility_exit_with_message("Can't load PDB rotamers from an empty filename.");
	}
	if ( !restype.is_ligand() ) { // Historical tracer output
		TR.Debug << "Warning: using PDB_ROTAMERS for non-ligand ResidueType!" << std::endl;
	}
	TR.Debug << "Initializing conformer library for " << rotamers_file << std::endl;
	SingleLigandRotamerLibraryOP pdb_rotamers( new SingleLigandRotamerLibrary() );
	pdb_rotamers->init_from_file( rotamers_file, restype );
	return pdb_rotamers;
}

std::string
SingleLigandRotamerLibraryCreator::keyname() const {
	return core::chemical::rotamers::PDBRotamerLibrarySpecification::library_name();
}

} //namespace rotamers
} //namespace pack
} //namespace core

