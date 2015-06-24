// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibraryCreator.cc
/// @brief  Class for instantiating a particular SingleResidueRotamerLibrary
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/rotamers/SingleResidueRotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Program header
#include <core/chemical/ResidueType.hh>

// Utility headers

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace rotamers {

/// @details Base class implementation ignores passed Residue object.
core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResidueRotamerLibraryCreator::create( core::chemical::ResidueType const & restype, core::conformation::Residue const & ) const {
	return create( restype );
}

} //namespace rotamers
} //namespace pack
} //namespace core

