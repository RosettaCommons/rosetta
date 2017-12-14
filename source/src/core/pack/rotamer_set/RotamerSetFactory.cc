// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSetFactory.hh
/// @brief  Residue Set Factory class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit header
#include <core/pack/rotamer_set/RotamerSetFactory.hh>

// Package headers
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
//#include <core/pack/rotamer_set/NucleicAcidRotamerSet.hh>

// Project headers
#include <core/conformation/Residue.hh>

// for symmetry
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/exit.hh>

// STL Headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

RotamerSetFactory::~RotamerSetFactory() = default;

RotamerSetOP
RotamerSetFactory::create_rotamer_set( conformation::Residue const & /*res*/ )
{

	// RM: Basing this on the commandline setting seems bad - what if we're setting symmetry from within RosettaScripts or the like?
	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() ) {
		return RotamerSetOP( new symmetry::SymmetricRotamerSet_() );
	}
	return RotamerSetOP( new RotamerSet_() );
}

}
}
}
