// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/dunbrack/SingleResidueRotamerLibrary.cc
/// @brief  SingleResidueRotamerLibrary class
/// @author Andrew Leaver-Fay
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit Headers
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>

#include <basic/Tracer.hh>

namespace core {
namespace pack {
namespace dunbrack {

static thread_local basic::Tracer TR( "core.pack.dunbrack.SingleResidueRotamerLibrary" );

SingleResidueRotamerLibrary::~SingleResidueRotamerLibrary()
{}

/// @brief Equality test for equivalence.
/// Two SingleResidueRotamerLibraries test equal if and only if they represent the exact same behavior
bool
SingleResidueRotamerLibrary::operator ==( SingleResidueRotamerLibrary const & ) const {
	// If you're comparing arbitrary SingleResidueRotamerLibrary, chances are they aren't equal.
	// (Override your subclass if this doesn't work for you.)
	TR.Warning << "[ WARNING ] Program is trying to compare two arbitrary SingleResidueRotamerLibraries - this is probably a bug." << std::endl;
	return false;
}

} // namespace dunbrack
} // namespace pack
} // namespace core
