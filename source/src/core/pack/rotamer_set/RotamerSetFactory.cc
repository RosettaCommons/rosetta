// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

RotamerSetFactory::~RotamerSetFactory() {}

RotamerSetOP
RotamerSetFactory::create_rotamer_set( conformation::Residue const & res )
{

	if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() ) {
		if ( res.is_protein() ) { // This check will be removed when we get rotamers for NAs and Ligands online
			return RotamerSetOP( new symmetry::SymmetricRotamerSet_() );
		} else {
			//std::cout << "[ WARNING ] PB HACK -- SHOULD DIE HERE!" << std::endl; // seems OK?
			return RotamerSetOP( new symmetry::SymmetricRotamerSet_() );
			//utility_exit_with_message( "Error in RotamerSetFactory, unsupported packing object" ); // get backtrace in gdb
			//exit(1); // add grace
			//return new AminoAcidRotamerSet(); // appease compiler
		}
	}

	if ( res.is_protein() ) { // This check will be removed when we get rotamers for NAs and Ligands online
		return RotamerSetOP( new RotamerSet_() );
	} else {
		//std::cout << "[ WARNING ] PB HACK -- SHOULD DIE HERE!" << std::endl; // seems OK?
		return RotamerSetOP( new RotamerSet_() );
		//utility_exit_with_message( "Error in RotamerSetFactory, unsupported packing object" ); // get backtrace in gdb
		//exit(1); // add grace
		//return new AminoAcidRotamerSet(); // appease compiler
	}
}

}
}
}
