// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/util.cc
/// @brief  Some useful nonmember functions

// Unit headers
#include <core/scoring/rna/util.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer tr( "core.scoring.rna.util" );

namespace core {
namespace scoring {
namespace rna {

Size
rna_residue_name_to_num( char const c )
{
	if ( c == 'a' ) return 1;
	if ( c == 'c' ) return 2;
	if ( c == 'g' ) return 3;
	if ( c == 'u' ) return 4;
	if ( c == 't' ) return 4;
	if ( c == 'Z' ) return 5; // Mg(2+)
	tr << "What is this? " << c << std::endl;
	utility_exit_with_message( "Asked for rna_residue_name_to_num for unknown residue_name" );
	return 0;
}

Size
protein_atom_name_to_num( std::string const & name )
{

	std::string blank_resname = "";
	return protein_atom_name_to_num( name, blank_resname );
}

Size
protein_atom_name_to_num( std::string const & name, std::string const & resname )
{
	// I should do this in a smarter way later, for now I'm just copying what is done for RNA

	if ( name == "Nbb" ) return 1;
	if ( name == " N  " ) return 1;
	if ( name == "CAbb" ) return 2;
	if ( name == " CA " ) return 2;
	if ( name == "CB" ) return 3;
	if ( name == " CB " ) return 3;
	if ( name == "CObb" ) return 4;
	if ( name == " C  " ) return 4;
	if ( name == "OCbb" ) return 5;
	if ( name == " O  " ) return 5;
	if ( name == " CEN" && resname != "" ) {
		if ( resname == "ALA" ) return 6;
		if ( resname == "CYS" ) return 7;
		if ( resname == "ASP" ) return 8;
		if ( resname == "GLU" ) return 9;
		if ( resname == "PHE" ) return 10;
		if ( resname == "GLY" ) return 11;
		if ( resname == "HIS" ) return 12;
		if ( resname == "ILE" ) return 13;
		if ( resname == "LYS" ) return 14;
		if ( resname == "LEU" ) return 15;
		if ( resname == "MET" ) return 16;
		if ( resname == "ASN" ) return 17;
		if ( resname == "PRO" ) return 18;
		if ( resname == "GLN" ) return 19;
		if ( resname == "ARG" ) return 20;
		if ( resname == "SER" ) return 21;
		if ( resname == "THR" ) return 22;
		if ( resname == "VAL" ) return 23;
		if ( resname == "TRP" ) return 24;
		if ( resname == "TYR" ) return 25;
	} else if ( resname == "" ) {
		if ( name == "CEN_ALA" ) return 6;
		if ( name == "CEN_CYS" ) return 7;
		if ( name == "CEN_ASP" ) return 8;
		if ( name == "CEN_GLU" ) return 9;
		if ( name == "CEN_PHE" ) return 10;
		if ( name == "CEN_GLY" ) return 11;
		if ( name == "CEN_HIS" ) return 12;
		if ( name == "CEN_ILE" ) return 13;
		if ( name == "CEN_LYS" ) return 14;
		if ( name == "CEN_LEU" ) return 15;
		if ( name == "CEN_MET" ) return 16;
		if ( name == "CEN_ASN" ) return 17;
		if ( name == "CEN_PRO" ) return 18;
		if ( name == "CEN_GLN" ) return 19;
		if ( name == "CEN_ARG" ) return 20;
		if ( name == "CEN_SER" ) return 21;
		if ( name == "CEN_THR" ) return 22;
		if ( name == "CEN_VAL" ) return 23;
		if ( name == "CEN_TRP" ) return 24;
		if ( name == "CEN_TYR" ) return 25;
	}

	std::cout << "What is this? " << name << std::endl;
	utility_exit_with_message( "Asked for protein_atom_name_to_num for unknown atom_name" );
	return 0;
}

} //rna
} //scoring
} //core
