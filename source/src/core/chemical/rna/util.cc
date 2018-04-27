// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rna/RNA_Util.cc
/// @author Rhiju Das

// Unit headers
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/types.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/Stub.hh>

// Project headers
#include <numeric/constants.hh>
#include <numeric/angle.functions.hh>

#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <ObjexxFCL/string.functions.hh>

// Utility headers

// C++

namespace core {
namespace chemical {
namespace rna {

static basic::Tracer TR( "core.chemical.rna" );

///////////////////////////////////////////////////////////////////////////////
Size
convert_acgu_to_1234( char const c )
{
	if ( rna_nts.find( c ) != std::string::npos ) {
		return rna_nts.find( c ) + 1;
	}
	if ( c == 't' ) return rna_nts.find( 'u' ) + 1;
	return 0;
}

/////////////////////////////////////////////////////////////////////////////
char get_edge_from_num( Size const num ) {
	if ( num == WATSON_CRICK ) return 'W';
	if ( num == HOOGSTEEN )    return 'H';
	if ( num == SUGAR )        return 'S';
	if ( num == O2PRIME )      return '2';
	if ( num == PHOSPHATE )    return 'P';
	return 'X';
}

/////////////////////////////////////////////////////////////////////////////
std::string  //Parin March 7, 2011
get_full_edge_from_num( Size const num ) {
	if ( num == WATSON_CRICK ) return "WC";
	if ( num == HOOGSTEEN )    return "HOOG";
	if ( num == SUGAR )        return "SUGAR";
	if ( num == O2PRIME )      return "O2PRIME";
	if ( num == PHOSPHATE )    return "PHOS";

	std::cout << "Invalid edge num = " << num << std::endl;
	utility_exit_with_message( "Invalid edge num!" );
	return "ERROR";
}

/////////////////////////////////////////////////////////////////////////////
BaseEdge
get_edge_from_char( char const e ) {
	if ( e == 'W' ) return WATSON_CRICK;
	if ( e == 'H' ) return HOOGSTEEN;
	if ( e == 'S' ) return SUGAR;
	if ( e == '2' ) return O2PRIME;
	if ( e == 'P' ) return PHOSPHATE;
	runtime_assert( e == 'X' );
	return ANY_BASE_EDGE;
}

/////////////////////////////////////////////////////////////////////////////
//This may be used elsewhere -- set up a util.hh?
char get_orientation_from_num( Size const num ) {
	if ( num == 1 ) return 'A';
	if ( num == 2 ) return 'P';
	return 'X';
}

/////////////////////////////////////////////////////////////////////////////
std::string
get_full_orientation_from_num( Size const num ) {
	if ( num == ANY_BASE_DOUBLET_ORIENTATION ) return "ANY ";
	if ( num == ANTIPARALLEL ) return "ANTI";
	if ( num == PARALLEL ) return "PARA";

	std::cout << "Invalid orientation num = " << num << std::endl;
	utility_exit_with_message( "Invalid orientation num!" );
	return "ERROR";
}

/////////////////////////////////////////////////////////////////////////////
char
get_side_from_num( Size const num ) {
	if ( num == ANY_BASE_STACK_SIDE ) return 'X';
	if ( num == ABOVE ) return 'A';
	if ( num == BELOW ) return 'B';
	utility_exit_with_message( "Invalid side num!" );
	return 'X';
}

/////////////////////////////////////////////////////////////////////////////
/// @details for stacking
std::string
get_full_side_from_num( Size const num ) {
	if ( num == ANY_BASE_STACK_SIDE ) return "ANY  ";
	if ( num == ABOVE ) return "ABOVE";
	if ( num == BELOW ) return "BELOW";
	utility_exit_with_message( "Invalid side num!" );
	return "ERROR";
}

/////////////////////////////////////////////////////////////////////////////
BaseDoubletOrientation
get_orientation_from_char( char const o ) {
	if ( o == 'A' ) return ANTIPARALLEL;
	if ( o == 'P' ) return PARALLEL;
	runtime_assert( o == 'X' );
	return ANY_BASE_DOUBLET_ORIENTATION;
}

/////////////////////////////////////////////////////////////////////////////
LW_BaseDoubletOrientation
get_LW_orientation_from_char( char const o ) {
	if ( o == 'C' ) return CIS;
	if ( o == 'T' ) return TRANS;
	runtime_assert( o == 'X' );
	return ANY_LW_BASE_DOUBLET_ORIENTATION;
}

/////////////////////////////////////////////////////////////////////////////
std::string //Parin April 19, 2011
get_full_LW_orientation_from_num( Size const num ){
	if ( num == 0 ) return "ANY  ";
	if ( num == 1 ) return "CIS  ";
	if ( num == 2 ) return "TRANS";

	std::cout << "Invalid orientation num = " << num << std::endl;
	utility_exit_with_message( "Invalid orientation num!" );
	return "ERROR";
}

/////////////////////////////////////////////////////////////////////////////
char
get_LW_orientation_from_num( Size const num ){
	if ( num == 0 ) return 'X';
	if ( num == 1 ) return 'C';
	if ( num == 2 ) return 'T';

	std::cout << "Invalid LW_orientation num = " << num << std::endl;
	utility_exit_with_message( "Invalid LW_orientation num!" );
	return 'X';
}

///////////////////////////////////////////////////////////////////////////////
std::string const first_base_atom( chemical::ResidueType const & rsd ) {
	// if (rsd.name1() == 'a' || rsd.name1() == 'g' )  return " N9 ";
	// return " N1 ";
	return rsd.atom_name( first_base_atom_index( rsd ) );
}


///////////////////////////////////////////////////////////////////////////////
Size first_base_atom_index( chemical::ResidueType const & rsd ) {
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 1 /*chi # 1 must be nucleic acid "chi"*/ );
	return atom_indices[ 3 ]; /* C2' ... C1' ... first base atom ...  chi1 torsion atom*/
}

///////////////////////////////////////////////////////////////////////////////
std::string const chi1_torsion_atom( chemical::ResidueType const & rsd ) {
	// if (rsd.name1() == 'a' || rsd.name1() == 'g' )  return " N9 ";
	// return " N1 ";
	return rsd.atom_name( chi1_torsion_atom_index( rsd ) );
}


///////////////////////////////////////////////////////////////////////////////
Size chi1_torsion_atom_index( chemical::ResidueType const & rsd ) {
	// HEY MAKE THIS MORE GENERAL? Maybe look at chi1?
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 1 /*chi # 1 must be nucleic acid "chi"*/ );
	return atom_indices[ 4 ]; /* C2' ... C1' ... first base atom ... chi1 torsion atom*/
}


///////////////////////////////////////////////////////////////////////////////
// consider moving this to chemical/util.cc.
std::string const default_jump_atom( chemical::ResidueType const & rsd ) {
	if ( rsd.is_RNA() || rsd.is_TNA() ) {
		if ( !rsd.is_coarse() ) {
			return chi1_torsion_atom( rsd );
		} else {
			return " Y  ";
		}
	} else if ( rsd.is_DNA() ) {
		return chi1_torsion_atom( rsd );
	}
	if ( rsd.is_carbohydrate() ) {
		return " C1 ";
	}
	if ( rsd.is_protein() ) {
		return " CA "; // note that this does not match 'traditional' choice in FoldTree.cc
	}
	if ( rsd.name3() == " MG" )  return "MG  ";
	if ( rsd.name3() == "HOH" )  return " O  ";
	if ( rsd.name3() == " ZN" )  return "ZN  ";
	if ( rsd.name3() == "XXX" ) return " Y  ";

	if ( rsd.name1() == 'Z' && !rsd.is_polymer() ) return rsd.atom_name( 1 ); // some other ligand

	std::cerr << "Residue ??? " << rsd.name3() << std::endl;
	utility_exit_with_message( "Do not know jump atom for this residue" );

	return "????";
}

///////////////////////////////////////////////////////////////////////////////
bool
possibly_canonical( chemical::AA const & aa1,  chemical::AA const & aa2 ) {
	using namespace core::chemical;
	return ( ( aa1 == na_rgu && aa2 == na_rcy ) ||
		( aa1 == na_rcy && aa2 == na_rgu ) ||
		( aa1 == na_rgu && aa2 == na_ura ) ||
		( aa1 == na_ura && aa2 == na_rgu ) ||
		( aa1 == na_rad && aa2 == na_ura ) ||
		( aa1 == na_ura && aa2 == na_rad )  );
}

///////////////////////////////////////////////////////////////////////////////
//no G-U
bool
possibly_canonical_strict( chemical::AA const & aa1,  chemical::AA const & aa2 ) {
	using namespace core::chemical;
	return ( ( aa1 == na_rgu && aa2 == na_rcy ) ||
		( aa1 == na_rcy && aa2 == na_rgu ) ||
		( aa1 == na_rad && aa2 == na_ura ) ||
		( aa1 == na_ura && aa2 == na_rad )  );
}

///////////////////////////////////////////////////////////////////////////
std::string get_WC_atom( core::chemical::AA const & res_type ){
	using namespace core::chemical;
	std::string WC_atom( "" );
	if ( res_type == na_rad ) WC_atom = " N1 ";
	if ( res_type == na_rcy ) WC_atom = " N3 ";
	if ( res_type == na_rgu ) WC_atom = " N1 ";
	if ( res_type == na_ura ) WC_atom = " N3 ";
	return WC_atom;
}

///////////////////////////////////////////////////////////////////////////////
void
get_watson_crick_base_pair_atoms(
	chemical::ResidueType const & rsd_type1,
	chemical::ResidueType const & rsd_type2,
	std::string & atom1,
	std::string & atom2 ) {

	using namespace core::chemical;

	// Defaults automatically to returning aa()
	AA const & aa1 = rsd_type1.na_analogue();
	AA const & aa2 = rsd_type2.na_analogue();

	if ( aa1 == na_rad && aa2 == na_ura ) {
		atom1 = " N1 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1 == na_ura && aa2 == na_rad ) {
		atom1 = " N3 ";
		atom2 = " N1 ";
		return;
	}
	if ( aa1 == na_rgu && aa2 == na_rcy ) {
		atom1 = " N1 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1 == na_rcy && aa2 == na_rgu ) {
		atom1 = " N3 ";
		atom2 = " N1 ";
		return;
	}
	if ( aa1 == na_rgu && aa2 == na_ura ) {
		atom1 = " O6 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1 == na_ura && aa2 == na_rgu ) {
		atom1 = " N3 ";
		atom2 = " O6 ";
		return;
	}

	atom1 = "XXXX";
	atom2 = "XXXX";
}

/////////////////////////////////////////////////////////////////////
void
get_base_pair_atoms(
	chemical::ResidueType const & rsd_type1,
	chemical::ResidueType const & rsd_type2,
	utility::vector1< std::string > & atom_ids1,
	utility::vector1< std::string > & atom_ids2,
	chemical::rna::BaseEdge const edge1,
	chemical::rna::BaseEdge const edge2,
	chemical::rna::LW_BaseDoubletOrientation const orientation
) {
	// AMW TODO: Cis BPs between the same AA are problematic.
	// (They're order dependent.)
	// AMW TODO: In 2018, convert this to some kind of map.

	// All possible BPs -- counting orientation -- so far observed in NDB:
	// cWW: All but GG; AA, CC, CU are rare
	// tWW: AA, GC, UA most common; AC, CC, UC, UG, UU, GG possible; AG unseen
	// cWH: AA, AC, CA, GC, GU, UC impossible; GG, UA most common; CC, GA, UG, UU, AG, AU, CU, CG possible
	// tWH: UA most common, then (still fairly common) AA, CA, GG; then CC, AG, CG, UG, GU, UU
	// cWS: all possible. AA, AC most common
	// tWS: AG, GU most common; only GA and GG are impossible
	// cHS: UG best, followed by CA, AA; all but GC, GU possible
	// tHS: AG best, then AA; all but GA, GC, GU, CG, UC, UU possible
	// cSS: all possible but AA, CA, UA, GA best
	// tSS: CC, UU impossible; AG most common, then AA, AC, AU

	using namespace chemical;

	atom_ids1.clear();
	atom_ids2.clear();

	// Enable na_analogue once it is also kept track of by BasePair
	//AA const aa1 = rsd_type1.aa();// == aa_unk ? rsd_type1.na_analogue() : rsd_type1.aa();
	//AA const aa2 = rsd_type2.aa();// == aa_unk ? rsd_type2.na_analogue() : rsd_type2.aa();

	// Enable base_analogue NOW and use PNA as an example
	AA const aa1 = ( rsd_type1.aa() == aa_unk || rsd_type1.aa() == aa_unp ) ? rsd_type1.base_analogue() : rsd_type1.aa();
	AA const aa2 = ( rsd_type2.aa() == aa_unk || rsd_type2.aa() == aa_unp ) ? rsd_type2.base_analogue() : rsd_type2.aa();

	if ( edge1 == WATSON_CRICK && edge2 == WATSON_CRICK ) {
		if ( orientation == CIS ) {
			if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H3 " );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O4 " );
			} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2 " );
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H41" );
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H3 " );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " O2 " );
			} else if ( aa2 == na_rad && aa1 == na_ura ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H3 " );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O4 " );
			} else if ( aa2 == na_rgu && aa1 == na_rcy ) {
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N3 " );
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2 " );
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H41" );
			} else if ( aa2 == na_rgu && aa1 == na_ura ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H3 " );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " O2 " );

				// special case -- isoG/isoC
			} else if ( rsd_type2.name3() == " IG" && rsd_type1.name3() == " IC" ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O4 " );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N3 " );
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H21" );
			} else if ( rsd_type1.name3() == " IG" && rsd_type2.name3() == " IC" ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O4 " );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H21" );


				// Noncanonical cWW base pairs
				// cWW: All but GG; AA, CC, CU are rare
			} else if ( aa1 == na_rad && aa2 == na_rad ) {
				// Could go either way?
				// BasePair needs a way to say "for this type of ,
				// I need a way to flip the atom IDs, but not for others (where that
				// wouldn't make sense because different edges or AAs are used."
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H61" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " N1 " ); // ugh
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " O2 " ); // ugh
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_rcy ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " O2 " ); // ugh
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H1 " );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O6 " );
			} else if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H1 " );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O6 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				// Could go either way?
				// BasePair needs a way to say "for this type of ,
				// I need a way to flip the atom IDs, but not for others (where that
				// wouldn't make sense because different edges or AAs are used."
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H41" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " N3 " ); // ugh
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O4 " );
			} else if ( aa2 == na_rcy && aa1 == na_ura ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O4 " );
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				// Could go either way?
				// BasePair needs a way to say "for this type of ,
				// I need a way to flip the atom IDs, but not for others (where that
				// wouldn't make sense because different edges or AAs are used."
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H3 " );
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " O2 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else { // trans
			// tWW: AA, GC, UA most common; AC, CC, UC, UG, UU, GG possible; AG unseen
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H61" );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N1 " );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H41" );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_rcy ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H41" );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H3 " );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O2 " );
			} else if ( aa2 == na_rad && aa1 == na_ura ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H3 " );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O2 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H41" );
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O2 " );
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H21" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H1 " );
			} else if ( aa2 == na_rcy && aa1 == na_rgu ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H21" );
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H1 " );
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O2 " );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H3 " );
			} else if ( aa2 == na_rcy && aa1 == na_ura ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O2 " );
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H3 " );
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " O2 " );
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H1 " );
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H3 " );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " O4 " );
			} else if ( aa2 == na_rgu && aa1 == na_ura ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H3 " );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " O4 " );
			} else if ( aa2 == na_ura && aa1 == na_ura ) {
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " O4 " );
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H3 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == WATSON_CRICK && edge2 == HOOGSTEEN ) {
		// cWH: AA, AC, CA, GC, GU, UC impossible; GG, UA most common; CC, GA, UG, UU, AG, AU, CU, CG possible
		// tWH: UA most common, then (still fairly common) AA, CA, GG; then CC, AG, CG, UG, GU, UU
		if ( orientation == CIS ) {
			if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O6 " );
				//atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " N7 " ); // Not confident in eventual length.
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H3 " );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O4 " );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) { // Part of the jump library
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H42" );
				// Not a good interaction but necessary to anchor this BP vs. trans
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " H41" );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H42" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H5 " ); // lame
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O6 " );
				//atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " N7 " ); // not really
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O4 " );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H5 " ); // lame
			} else if ( aa1 == na_rgu && aa2 == na_rad ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H62" );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " N7 " );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " O6 " );
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H62" );
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H8 " ); // lame
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H5 " ); // lame
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " O4 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H62" );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N7 " );
				//atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " O6 " ); // can't trust distance
			} else if ( aa1 == na_rad && aa2 == na_rcy ) { // Part of the jump library
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H42" );
				// Not a good interaction but necessary to anchor this BP vs. trans
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " H5 " );
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H62" );
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H5 " );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H42" );
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				//atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " 06 " ); // can't trust
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N7 " );
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O6 " );
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				atom_ids1.push_back( " 06 " );  atom_ids2.push_back( " H5 " );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " 04" );
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " N7 " );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H62" );
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H8 " ); // lame
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H5 " ); // lame
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " O4 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == HOOGSTEEN && edge2 == WATSON_CRICK ) {
		// cWH: AA, AC, CA, GC, GU, UC impossible; GG, UA most common; CC, GA, UG, UU, AG, AU, CU, CG possible
		// tWH: UA most common, then (still fairly common) AA, CA, GG; then CC, AG, CG, UG, GU, UU
		if ( orientation == CIS ) {
			if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O6 " );
				//atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " N7 " ); // Not confident in eventual length.
			} else if ( aa2 == na_rad && aa1 == na_ura ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H3 " );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O4 " );
			} else if ( aa2 == na_rad && aa1 == na_rcy ) { // Part of the jump library
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H42" );
				// Not a good interaction but necessary to anchor this BP vs. trans
				atom_ids2.push_back( " H2 " );  atom_ids1.push_back( " H41" );
			} else if ( aa2 == na_rcy && aa1 == na_rcy ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H42" );
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H5 " ); // lame
			} else if ( aa2 == na_rcy && aa1 == na_rgu ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O6 " );
				//atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " N7 " ); // not really
			} else if ( aa2 == na_rcy && aa1 == na_ura ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O4 " );
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H5 " ); // lame
			} else if ( aa2 == na_rgu && aa1 == na_rad ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H62" );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_rgu && aa1 == na_rgu ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " N7 " );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " O6 " );
			} else if ( aa2 == na_ura && aa1 == na_rad ) {
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H62" );
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_ura && aa1 == na_rgu ) {
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H8 " ); // lame
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_ura && aa1 == na_ura ) {
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H5 " ); // lame
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " O4 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa2 == na_rad && aa1 == na_rad ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H62" );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N7 " );
				//atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " O6 " ); // can't trust distance
			} else if ( aa1 == na_rad && aa2 == na_rcy ) { // Part of the jump library
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H42" );
				// Not a good interaction but necessary to anchor this BP vs. trans
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " H5 " );
			} else if ( aa2 == na_rcy && aa1 == na_rad ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H62" );
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_rcy && aa1 == na_rcy ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H5 " );
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H42" );
			} else if ( aa2 == na_rcy && aa1 == na_rgu ) {
				//atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " 06 " ); // can't trust
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_rgu && aa1 == na_rgu ) {
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N7 " );
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O6 " );
			} else if ( aa2 == na_rgu && aa1 == na_ura ) {
				atom_ids2.push_back( " 06 " );  atom_ids1.push_back( " H5 " );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " 04" );
			} else if ( aa2 == na_ura && aa1 == na_rad ) {
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " N7 " );
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H62" );
			} else if ( aa2 == na_ura && aa1 == na_rgu ) {
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( " H8 " ); // lame
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " N7 " );
			} else if ( aa2 == na_ura && aa1 == na_ura ) {
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( " H5 " ); // lame
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " O4 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == HOOGSTEEN && edge2 == HOOGSTEEN ) {
		// cHH: GA, GC, GG: possible but rare
		// tHH: all possible but CC, UU, GU; AA by far most common
		if ( orientation == CIS ) {
			if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O6 " );
			} else if ( aa1 == na_rgu && aa2 == na_rad ) {
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O6 " );
			} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H6 " ); // lame
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H6 " ); // lame
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				// Asymmetric, of course.
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H8 " ); // lame
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}

		} else {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " N7 " );  atom_ids2.push_back( " H62" );
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( " N7 " );  atom_ids2.push_back( " H42" );
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				atom_ids2.push_back( " N7 " );  atom_ids1.push_back( " H42" );
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O6 " );
			} else if ( aa1 == na_rgu && aa2 == na_rad ) {
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O6 " );
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O4 " );
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O4 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " N7 " );
			} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " N7 " );
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O4 " );
			} else if ( aa1 == na_ura && aa2 == na_rcy ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O4 " );
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				// heck this is asymmetric too
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H8 " );
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == WATSON_CRICK && edge2 == SUGAR ) {
		if ( orientation == CIS ) {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H61" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H61" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				// these (and below) are VASTLY different cis vs. trans.
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H41" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H41" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rgu && aa2 == na_rad ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
				//atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
				//atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " H6 " ); // nah
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H2 " ); // lame
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O2'" );
				//atom_ids1.push_back( " H61" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H22" );
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H61" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O2'" );
				//atom_ids1.push_back( " H41" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O2'" );
				//atom_ids1.push_back( " H41" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_ura && aa2 == na_rcy ) {
				//atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " H6 " ); // nah
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( "HO2'" );
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				//atom_ids1.push_back( " H3 " );  atom_ids2.push_back( " H6 " ); // nah
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == SUGAR && edge2 == WATSON_CRICK ) {
		if ( orientation == CIS ) {
			if ( aa2 == na_rad && aa1 == na_rad ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( "HO2'" );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_rcy ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H61" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( "HO2'" );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_ura ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H61" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rcy && aa1 == na_rad ) {
				// these (and below) are VASTLY different cis vs. trans.
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( "HO2'" );
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rcy && aa1 == na_rcy ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H41" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rcy && aa1 == na_rgu ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( "HO2'" );
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rcy && aa1 == na_ura ) {
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H41" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rgu && aa1 == na_rad ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N3 " );
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H2 " ); // lame
			} else if ( aa2 == na_rgu && aa1 == na_rcy ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2'" );
				//atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rgu && aa1 == na_rgu ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H22" );
			} else if ( aa2 == na_rgu && aa1 == na_ura ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2'" );
				//atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " H6 " ); // nah
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa2 == na_rad && aa1 == na_rad ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H2 " ); // lame
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_rcy ) {
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O2'" );
				//atom_ids2.push_back( " H61" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H22" );
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_ura ) {
				atom_ids2.push_back( " N1 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H61" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rcy && aa1 == na_rad ) {
				// Could require "O4 atom name" for thio residues.
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O4'" );
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rcy && aa1 == na_rcy ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O2'" );
				//atom_ids2.push_back( " H41" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rcy && aa1 == na_rgu ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rcy && aa1 == na_ura ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O2'" );
				//atom_ids2.push_back( " H41" );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rgu && aa1 == na_rcy ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_rgu && aa1 == na_ura ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( "HO2'" );
				//atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " H6 " ); // nah
			} else if ( aa2 == na_ura && aa1 == na_rad ) {
				atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " N3 " );
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H2 " ); // lame
			} else if ( aa2 == na_ura && aa1 == na_rcy ) {
				//atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " H6 " ); // nah
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( "HO2'" );
			} else if ( aa2 == na_ura && aa1 == na_rgu ) {
				atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H22" );
			} else if ( aa2 == na_ura && aa1 == na_ura ) {
				//atom_ids2.push_back( " H3 " );  atom_ids1.push_back( " H6 " ); // nah
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( "HO2'" ); // lame
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == HOOGSTEEN && edge2 == SUGAR ) {
		if ( orientation == CIS ) {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " N7 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) { // part of jump library
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O2 " );
				atom_ids1.push_back( " N7 " );  atom_ids2.push_back( "HO2'" );
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rgu && aa2 == na_rad ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_ura && aa2 == na_rcy ) {
				atom_ids1.push_back( " H5 " );  atom_ids2.push_back( " O2 " );
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H22" ); // lame
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				atom_ids1.push_back( " H5 " );  atom_ids2.push_back( " O2 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " N7 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O2 " );
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " N7 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H62" );  atom_ids2.push_back( " O2 " );
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O2 " );
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( " H41" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " H42" );  atom_ids2.push_back( " O2 " );
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				atom_ids1.push_back( " O4 " );  atom_ids2.push_back( " H22" );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == SUGAR && edge2 == HOOGSTEEN ) {
		if ( orientation == CIS ) {
			if ( aa2 == na_rad && aa1 == na_rad ) {
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " N7 " );  atom_ids1.push_back( " H22" );
			} else if ( aa2 == na_rad && aa1 == na_rcy ) { // part of jump library
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O2 " );
				atom_ids2.push_back( " N7 " );  atom_ids1.push_back( "HO2'" );
			} else if ( aa2 == na_rcy && aa1 == na_rad ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rcy && aa1 == na_rgu ) {
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rgu && aa1 == na_rad ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H2 " ); // lame
			} else if ( aa2 == na_rgu && aa1 == na_rgu ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H22" );
			} else if ( aa2 == na_ura && aa1 == na_rad ) {
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( " H2 " ); // lame
			} else if ( aa2 == na_ura && aa1 == na_rcy ) {
				atom_ids2.push_back( " H5 " );  atom_ids1.push_back( " O2 " );
			} else if ( aa2 == na_ura && aa1 == na_rgu ) {
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( " H22" ); // lame
			} else if ( aa2 == na_ura && aa1 == na_ura ) {
				atom_ids2.push_back( " H5 " );  atom_ids1.push_back( " O2 " );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa2 == na_rad && aa1 == na_rad ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " N3 " );
				atom_ids2.push_back( " N7 " );  atom_ids1.push_back( " H2 " ); // lame
			} else if ( aa2 == na_rad && aa1 == na_rcy ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O2 " );
			} else if ( aa2 == na_rad && aa1 == na_rgu ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " N3 " );
				atom_ids2.push_back( " N7 " );  atom_ids1.push_back( " H22" );
			} else if ( aa2 == na_rad && aa1 == na_ura ) {
				atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H62" );  atom_ids1.push_back( " O2 " );
			} else if ( aa2 == na_rcy && aa1 == na_rad ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " N3 " );
			} else if ( aa2 == na_rcy && aa1 == na_rcy ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O2 " );
			} else if ( aa2 == na_rcy && aa1 == na_ura ) {
				atom_ids2.push_back( " H41" );  atom_ids1.push_back( " O2'" );
				atom_ids2.push_back( " H42" );  atom_ids1.push_back( " O2 " );
			} else if ( aa2 == na_rgu && aa1 == na_rgu ) {
				atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H22" );
			} else if ( aa2 == na_ura && aa1 == na_rad ) {
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( " H2 " ); // lame
			} else if ( aa2 == na_ura && aa1 == na_rgu ) {
				atom_ids2.push_back( " O4 " );  atom_ids1.push_back( " H22" );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else if ( edge1 == SUGAR && edge2 == SUGAR ) {
		if ( orientation == CIS ) {
			// Obviously -- symmetry issues.
			if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " N3 " ); // lame
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " O2 " ); // lame
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " N3 " ); // lame
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " O2 " ); // lame
			} else if ( aa1 == na_rcy && aa2 == na_rad ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_rcy && aa2 == na_rcy ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_rcy && aa2 == na_rgu ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_rcy && aa2 == na_ura ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_rad && aa2 == na_rad ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H22" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_rcy ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H22" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H22" );  atom_ids2.push_back( " N3 " );
			} else if ( aa1 == na_rad && aa2 == na_ura ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( "HO2'" );
				//atom_ids1.push_back( " H22" );  atom_ids2.push_back( " H6 " ); // nah
			} else if ( aa1 == na_ura && aa2 == na_rad ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_ura && aa2 == na_rcy ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_ura && aa2 == na_rgu ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else if ( aa1 == na_ura && aa2 == na_ura ) {
				atom_ids1.push_back( "HO2'" );  atom_ids2.push_back( " O2'" );
				atom_ids1.push_back( " O2 " );  atom_ids2.push_back( "HO2'" ); // lame
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		} else {
			if ( aa1 == na_rad && aa2 == na_rad ) {
				// Some of these interactions must be transient. The N1-HO2' makes it
				// more like a WC-S, so I'm going for the reflexive one
				//atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " N3 " ); // lame
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H2 " ); // lame
			} else if ( aa1 == na_rad && aa2 == na_rgu ) {
				// Some of these interactions must be transient. The N1-HO2' makes it
				// more like a WC-S, so I'm going for the reflexive one
				//atom_ids1.push_back( " N1 " );  atom_ids2.push_back( "HO2'" );
				atom_ids1.push_back( " H2 " );  atom_ids2.push_back( " N3 " ); // lame
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H22" );
			} else if ( aa1 == na_rgu && aa2 == na_rad ) {
				// Some of these interactions must be transient. The N1-HO2' makes it
				// more like a WC-S, so I'm going for the reflexive one
				//atom_ids2.push_back( " N1 " );  atom_ids1.push_back( "HO2'" );
				atom_ids2.push_back( " H2 " );  atom_ids1.push_back( " N3 " ); // lame
				atom_ids2.push_back( " N3 " );  atom_ids1.push_back( " H22" );
			} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
			} else if ( aa2 == na_rgu && aa1 == na_rcy ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2'" );
			} else if ( aa1 == na_rgu && aa2 == na_rgu ) {
				atom_ids1.push_back( " O2'" );  atom_ids2.push_back( " H21" );
				atom_ids1.push_back( " N3 " );  atom_ids2.push_back( " H22" );
				atom_ids1.push_back( " H22" );  atom_ids2.push_back( " N3 " );
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
			} else if ( aa1 == na_rgu && aa2 == na_ura ) {
				atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2'" );
			} else if ( aa2 == na_rgu && aa1 == na_ura ) {
				atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2'" );
			} else {
				TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
				TR.Debug << "These residues don't cleanly base pair in that orientation." << std::endl;
			}
		}
	} else {
		TR.Debug << "Requested " << edge1 << " " << edge2 << " base pairing atoms for " << aa1 << " " << aa2 << "." << std::endl;
		TR.Debug << "These residues might base pair in that orientation, but we don't support it yet." << std::endl;

	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
//Unify the version in StepWiseRNA_Utill.cc and RNA_CentroidInfo.cc on June 25, 2011
// Also, this copies some code from Phil's dna/base_geometry.cc
//Comments (Parin Sep 23,2009)...possible problem if every atoms in the nucleotide is virtual...in that case numatoms=0....will this crash the code??

Vector
get_rna_base_centroid( conformation::Residue const & rsd, bool verbose ){

	// The centroid is C for pdb_GAI -- we sort of need this by
	// definition. So that we don't refer to pdb_GAI by name,
	// let's say that all non-polymer residues have a centroid of
	// their first atom's xyz. This is vaguely reasonable a priori,
	// since it's slightly more likely than not to be the root, maybe,
	// and may prompt a better general centroid function.

	// Really, AMW TODO don't call get_rna_base_centroid on non RNA.
	// instead call something else!
	if ( !rsd.is_polymer() ) {
		Vector centroid( 0.0 );
		for ( Size ii = 1; ii <= rsd.nheavyatoms(); ++ii ) {
			centroid += rsd.xyz( ii );
		}
		return centroid / rsd.nheavyatoms();
		//return rsd.xyz(1);
	}

	// Rats... what if it is a 'polymer residue' like protein but not
	// in a polymer context?
	if ( rsd.is_polymer() && !rsd.has_lower_connect() && !rsd.has_upper_connect() ) {
		return rsd.xyz(1);
	}

	//SML PHENIX conference
	if ( !rsd.is_NA() ) {
		std::cout << "name " << rsd.type().name() << std::endl;
		if ( rsd.is_carbohydrate() || basic::options::option[basic::options::OptionKeys::rna::erraser::rna_prot_erraser].value() ) {
			return Vector( 0.0, 0.0, 0.0 );
		} else { //if not option
			utility_exit_with_message( "non - RNA residue inside get_rna_base_centroid" );
		}
	} //if not RNA

	Vector centroid( 0.0 );
	Size numatoms = 0;
	//std::cout << "mid of get centroid: " << centroid.x() << " " << centroid.y() << " "  << centroid.z() << std::endl;

	//Consistency check:
	//if(rsd.type().atom_name(rsd.first_sidechain_atom()) !=" O2'") utility_exit_with_message( "rsd.type().atom_name(rsd.first_sidechain_atom()) !=\" O2'\" " );
	//if(rsd.atom_name( rsd.first_sidechain_atom() )!=" O2'") utility_exit_with_message("rsd.atom_name( rsd.first_sidechain_atom() )!=\" O2'\"");

	// AMW: added is_RNA() condition because DNA may get down here.
	if ( rsd.is_RNA() && rsd.RNA_info().o2prime_index() != rsd.first_sidechain_atom() ) {
		utility_exit_with_message( "rsd.RNA_info().o2prime_index() != rsd.first_sidechain_atom() for residue "+rsd.name() );
	}

	if ( verbose )  std::cout << "Base atoms" << std::endl;

	for ( Size i = rsd.first_sidechain_atom() + 1; i <= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.

		if ( verbose ) std::cout << "atom " << i  << " " <<  "name = " << rsd.type().atom_name( i ) << " type = " << rsd.atom_type( i ).name()  << " " << rsd.atom_type_index( i ) << " " << rsd.atomic_charge( i );

		if ( rsd.is_RNA() && rsd.RNA_info().atom_is_virtual( i ) ) {
			if ( verbose ) std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}

		if ( rsd.is_RNA() && rsd.atom_type( i ).is_repulsive() ) {
			if ( verbose ) std::cout << "  Repulsive: Ignore! " << std::endl;
			continue;
		}

		if ( verbose ) std::cout << std::endl;

		centroid += rsd.xyz( i );
		//std::cout << "rsd " << rsd.seqpos() << " ATOM " << rsd.atom_name( i ) << centroid.x() << " " << centroid.y() << " "  << centroid.z() << std::endl;

		numatoms++;
	}

	// kludge... where are these atoms going???
	if ( numatoms == 0 || centroid.x() > 100000 || centroid.x() < -100000 ) { //Centroid not well defined in this case...probably because rsd is a virtual residue...just return 0
		Vector dummy_centroid( 0.0 );
		return dummy_centroid;
	}

	centroid /= static_cast< Real > ( numatoms );
	//std::cout << "bottom of get centroid: " << centroid.x() << " " << centroid.y() << " "  << centroid.z() << std::endl;

	return centroid;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Unify the version in StepWiseRNA_Util.cc and RNA_CentroidInfo.cc on June 25, 2011
numeric::xyzMatrix< core::Real >
get_rna_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ){

	//std::cout << "in get_rna_base_coordinate_system head Residue number " << rsd.seqpos() << " name " << rsd.type().name() << std::endl;
	//std::cout << "centroid: " << centroid.x() << " " << centroid.y() << " "  << centroid.z() << std::endl;
	// AMW TODO: similarly, don't call this RNA FUNCTION for non RNA residues.
	// I suppose ligands with planar bits with H-bond donors/acceptors can have
	// 'edges' in the same way, perhaps. Separate function!

	// We have to retain our mentions of pdb_GAI here because we are specialized
	// to its particular WC atoms.

	using namespace chemical;

	//SML PHENIX conference
	if ( !rsd.is_NA() && rsd.type().name() != "pdb_GAI" ) {
		if ( rsd.is_carbohydrate() || basic::options::option[basic::options::OptionKeys::rna::erraser::rna_prot_erraser].value() ) {
			return numeric::xyzMatrix< core::Real > ::identity();
		} else if ( rsd.is_polymer() && !rsd.has_lower_connect() && !rsd.has_upper_connect() ) {
			return numeric::xyzMatrix< core::Real > ::identity();
		} else if ( rsd.is_ligand() ) {
			return numeric::xyzMatrix< core::Real > ::identity();
		} else { //if not option
			utility_exit_with_message( "non - RNA residue inside get_rna_base_coordinate_system, abort" );
		}
	} //if not RNA

	AA res_type = rsd.aa();
	if ( res_type == aa_unk || res_type == aa_unp ) {
		// Prioritize base_analogue
		if ( rsd.type().base_analogue() != aa_unk &&
				rsd.type().base_analogue() != aa_unp ) {
			res_type = rsd.type().base_analogue();
		} else {
			res_type = rsd.type().na_analogue();
		}
	}
	// convert DNA to RNA equivalent
	if ( res_type == na_ade ) res_type = na_rad;
	if ( res_type == na_cyt ) res_type = na_rcy;
	if ( res_type == na_gua ) res_type = na_rgu;
	if ( res_type == na_thy ) res_type = na_ura;

	Vector x, y, z;

	Vector WC_coord, H_coord;
	if ( res_type == aa_unk || res_type == aa_unp ) {
		if ( rsd.type().name() == "pdb_GAI" ) {
			WC_coord = rsd.xyz("N2");
			H_coord = rsd.xyz("N3");
		} else {
			// Just use the first two sidechain atoms for generality
			// AMW TODO: instead you could imagine taking the else clause and permitting for is_purine or is_pyrimidine
			// (if no na_analogue)
			WC_coord = rsd.xyz( rsd.first_sidechain_atom() );
			H_coord = rsd.xyz( rsd.first_sidechain_atom() + 1 );
		}
	} else {
		// Make an axis pointing from base centroid to Watson-Crick edge.
		std::string WC_atom;
		if ( res_type == na_rad || res_type == na_rgu || res_type == na_lra || res_type == na_lrg ) WC_atom = "N1";
		if ( res_type == na_rcy || res_type == na_ura || res_type == na_lrc || res_type == na_lur ) WC_atom = "N3";
		// Make a perpendicular axis pointing from centroid towards
		// Hoogstein edge (e.g., major groove in a double helix).
		std::string H_atom;
		if ( res_type == na_rad || res_type == na_rgu || res_type == na_lra || res_type == na_lrg ) H_atom = "N7" ;
		if ( !rsd.has( H_atom ) && ( rsd.name3() == "7DA" || rsd.name3() == "FA7" ) ) H_atom = "C7";
		if ( res_type == na_rcy || res_type == na_ura || res_type == na_lrc || res_type == na_lur ) H_atom = "C5";
		if ( res_type == na_ura && rsd.name3() == "PSU" ) H_atom = "N1"; // pseudoU is flipped
		//std::cout << "WC_atom " << WC_atom << std::endl;

		WC_coord = rsd.xyz( WC_atom );
		H_coord = rsd.xyz( H_atom );
	}

	//std::cout << "WC_coord: " << WC_coord.x() << " " << WC_coord.y() << " "  << WC_coord.z() << std::endl;

	x = WC_coord - centroid;
	//std::cout << "x: " << x.x() << " " << x.y() << " "  << x.z() << std::endl;
	x.normalize();

	y = H_coord - centroid; //not orthonormal yet...
	//std::cout << "y: " << y.x() << " " << y.y() << " "  << y.z() << std::endl;
	z = cross( x, y );
	//std::cout << "z: " << z.x() << " " << z.y() << " "  << z.z() << std::endl;
	z.normalize(); // Should point roughly 5' to 3' if in a double helix.

	y = cross( z, x );
	//std::cout << "y: " << y.x() << " " << y.y() << " "  << y.z() << std::endl;
	y.normalize(); //not necessary but doesn't hurt.

	numeric::xyzMatrix< core::Real > M = numeric::xyzMatrix< core::Real > ::cols( x, y, z );
	return M;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief wrapper around get_rna_base_centroid & get_rna_base_coordinate_system
core::kinematics::Stub
get_rna_base_coordinate_system_stub( core::conformation::Residue const & rsd )
{
	Vector const centroid = get_rna_base_centroid( rsd, false /*verbose*/ );
	numeric::xyzMatrix< Real > const M = get_rna_base_coordinate_system( rsd, centroid );
	return core::kinematics::Stub( M, centroid );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Check whether one of the atom beyond to a base and another to a phosphate. Doesn't have to be on the same nucleotide.
bool
is_base_phosphate_atom_pair( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const atomno_1, Size const atomno_2 ){

	if ( !rsd_1.is_RNA() ) return false;
	if ( !rsd_2.is_RNA() ) return false;

	bool is_base_phosphate_atom_pair = false;

	if ( ( rsd_1.RNA_info().atom_is_phosphate( atomno_1 ) && ( rsd_2.RNA_info().is_RNA_base_atom( atomno_2 ) ) ) ) is_base_phosphate_atom_pair = true;
	if ( ( rsd_2.RNA_info().atom_is_phosphate( atomno_2 ) && ( rsd_1.RNA_info().is_RNA_base_atom( atomno_1 ) ) ) ) is_base_phosphate_atom_pair = true;

	if ( is_base_phosphate_atom_pair ) { //This Assume that rsd_1 and rsd_2 are the same!!!
		if ( rsd_1.seqpos() == rsd_2.seqpos() && ( rsd_1.path_distance( atomno_1, atomno_2 ) < 4 ) ) { //consistency check!
			utility_exit_with_message( "is_base_phosphate_atom_pair but rsd.path_distance( " + ObjexxFCL::string_of( atomno_1 ) + ", " + ObjexxFCL::string_of( atomno_2 ) + " ) < 4" );
		}
	}

	return is_base_phosphate_atom_pair;

}


//////////////////////////////////////////////////////////////////////////////////////
ChiState
get_residue_base_state( conformation::Residue const & rsd ) {
	Real const CHI_CUTOFF = 15.0; //Kinda RANDOM..ROUGH average between north chi_anti (~79) and north chi_syn (~-50)
	Real const chi = numeric::principal_angle_degrees( rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );
	if ( chi <= CHI_CUTOFF ) {
		return SYN;
	} else {
		return ANTI;
	}
}

//////////////////////////////////////////////////////////////////////////////////////
PuckerState
get_residue_pucker_state( conformation::Residue const & rsd ) {
	static RNA_FittedTorsionInfo const rna_fitted_torsion_info;
	Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );
	Real delta = numeric::principal_angle_degrees( rsd.mainchain_torsion( DELTA ) );
	if ( delta <= DELTA_CUTOFF ) {
		return NORTH;
	} else {
		return SOUTH;
	}
}



bool
rna_dna_match( core::chemical::AA const & aa1, core::chemical::AA const & aa2 )
{
	if ( aa1 == aa2 ) return true;
	using namespace core::chemical;
	switch ( aa1 ) {
	case( na_rad ) :
		return ( aa2 == na_ade );
	case( na_rcy ) :
		return ( aa2 == na_cyt );
	case( na_rgu ) :
		return ( aa2 == na_gua );
	case( na_ura ) :
		return ( aa2 == na_thy );
	case( na_ade ) :
		return ( aa2 == na_rad );
	case( na_cyt ) :
		return ( aa2 == na_rcy );
	case( na_gua ) :
		return ( aa2 == na_rgu );
	case( na_thy ) :
		return ( aa2 == na_ura );
	default :
		break;
	}
	return false;
}

} //ns rna
} //ns chemical
} //ns core
