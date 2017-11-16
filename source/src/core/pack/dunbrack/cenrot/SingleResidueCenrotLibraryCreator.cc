// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/cenrot/SingleResidueCenrotLibraryCreator.cc
/// @brief  Class for instantiating a particular SingleResidueRotamerLibrary
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibraryCreator.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Program header
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/CenrotRotamerLibrarySpecification.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
#include <core/pack/dunbrack/cenrot/CenrotLibrary.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

static basic::Tracer TR("core.pack.dunbrack.cenrot.SingleResidueCenrotLibraryCreator");

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueCenrotLibraryCreator::create( core::chemical::ResidueType const & restype) const {
	using namespace core::chemical::rotamers;
	RotamerLibrarySpecificationCOP libspec( restype.rotamer_library_specification() );
	// If the factory system is sound, these two checks should work.
	debug_assert( libspec );
	CenrotRotamerLibrarySpecificationCOP cenrot_libspec = utility::pointer::dynamic_pointer_cast< CenrotRotamerLibrarySpecification const >(libspec);
	debug_assert( cenrot_libspec );

	core::chemical::AA aa( cenrot_libspec->get_aa() );
	debug_assert( aa <= core::chemical::num_canonical_aas );
	core::pack::rotamers::SingleResidueRotamerLibraryCOP lib( core::pack::dunbrack::cenrot::CenrotLibrary::get_instance()->get_cenrot_library_by_aa(aa) );
	return lib;
}

std::string
SingleResidueCenrotLibraryCreator::keyname() const {
	return core::chemical::rotamers::CenrotRotamerLibrarySpecification::library_name();
}

} //namespace cenrot
} //namespace rotamers
} //namespace pack
} //namespace core

