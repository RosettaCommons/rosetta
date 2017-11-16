// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/dunbrack/SingleResidueDunbrackLibraryCreator.hh
/// @brief  Class for instantiating a particular SingleResidueRotamerLibrary
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Program header
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/DunbrackRotamerLibrarySpecification.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace dunbrack {

static basic::Tracer TR("core.pack.dunbrack.SingleResidueDunbrackLibraryCreator");

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueDunbrackLibraryCreator::create( core::chemical::ResidueType const & restype) const {
	using namespace core::chemical::rotamers;
	RotamerLibrarySpecificationCOP libspec( restype.rotamer_library_specification() );
	// If the factory system is sound, these two checks should work.
	debug_assert( libspec );
	DunbrackRotamerLibrarySpecificationCOP dun_libspec = utility::pointer::dynamic_pointer_cast< DunbrackRotamerLibrarySpecification const >(libspec);
	debug_assert( dun_libspec );

	core::chemical::AA aa( dun_libspec->get_aa() );
	// TODO: Split D_aas into their own Rotamer library set (wrapper of l-aas)
	if ( core::chemical::is_canonical_D_aa( aa ) ) {
		aa = core::chemical::get_L_equivalent( aa );
	}
	debug_assert( aa <= core::chemical::num_canonical_aas );
	core::pack::rotamers::SingleResidueRotamerLibraryCOP lib( core::pack::dunbrack::RotamerLibrary::get_instance()->get_library_by_aa(aa) );
	return lib;
}

std::string
SingleResidueDunbrackLibraryCreator::keyname() const {
	return core::chemical::rotamers::DunbrackRotamerLibrarySpecification::library_name();
}

} //namespace dunbrack
} //namespace pack
} //namespace core

