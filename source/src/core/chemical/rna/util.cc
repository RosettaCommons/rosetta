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
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <ObjexxFCL/string.functions.hh>

// Utility headers

// C++

namespace core {
namespace chemical {
namespace rna {

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
	if ( rsd.is_RNA() ) {
		if ( !rsd.is_coarse() ) {
			return chi1_torsion_atom( rsd );
		} else {
			return " Y  ";
		}
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
get_watson_crick_base_pair_atoms(
	chemical::ResidueType const & rsd_type1,
	chemical::ResidueType const & rsd_type2,
	utility::vector1< std::string > & atom_ids1,
	utility::vector1< std::string > & atom_ids2  )
{
	using namespace chemical;

	atom_ids1.clear();
	atom_ids2.clear();

	AA const & aa1 = rsd_type1.aa();
	AA const & aa2 = rsd_type2.aa();

	if ( aa1 == na_rad && aa2 == na_ura ) {
		atom_ids1.push_back( " N1 " );  atom_ids2.push_back( " H3 " );
		atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O4 " );
		return;
	} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
		atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N3 " );
		atom_ids1.push_back( " H21" );  atom_ids2.push_back( " O2 " );
		atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H41" );
		return;
	} else if ( aa1 == na_rgu && aa2 == na_ura ) {
		atom_ids1.push_back( " O6 " );  atom_ids2.push_back( " H3 " );
		atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " O2 " );
		return;
	} else if ( aa2 == na_rad && aa1 == na_ura ) {
		atom_ids2.push_back( " N1 " );  atom_ids1.push_back( " H3 " );
		atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O4 " );
		return;
	} else if ( aa2 == na_rgu && aa1 == na_rcy ) {
		atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N3 " );
		atom_ids2.push_back( " H21" );  atom_ids1.push_back( " O2 " );
		atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H41" );
		return;
	} else if ( aa2 == na_rgu && aa1 == na_ura ) {
		atom_ids2.push_back( " O6 " );  atom_ids1.push_back( " H3 " );
		atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " O2 " );
		return;

		// special case -- isoG/isoC
	} else if ( rsd_type2.name3() == " IG" && rsd_type1.name3() == " IC" ) {
		atom_ids2.push_back( " H61" );  atom_ids1.push_back( " O4 " );
		atom_ids2.push_back( " H1 " );  atom_ids1.push_back( " N3 " );
		atom_ids2.push_back( " O2 " );  atom_ids1.push_back( " H21" );
		return;
	} else if ( rsd_type1.name3() == " IG" && rsd_type2.name3() == " IC" ) {
		atom_ids1.push_back( " H61" );  atom_ids2.push_back( " O4 " );
		atom_ids1.push_back( " H1 " );  atom_ids2.push_back( " N3 " );
		atom_ids1.push_back( " O2 " );  atom_ids2.push_back( " H21" );
		return;
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////
//Unify the version in StepWiseRNA_Utill.cc and RNA_CentroidInfo.cc on June 25, 2011
// Also, this copies some code from Phil's dna/base_geometry.cc
//Comments (Parin Sep 23,2009)...possible problem if every atoms in the nucleotide is virtual...in that case numatoms=0....will this crash the code??

Vector
get_rna_base_centroid( conformation::Residue const & rsd, bool verbose ){

	//SML PHENIX conference
	if ( !rsd.is_RNA() ) {
		if ( basic::options::option[basic::options::OptionKeys::rna::erraser::rna_prot_erraser].value() ) {
			return Vector( 0.0, 0.0, 0.0 );
		} else { //if not option
			utility_exit_with_message( "non - RNA residue inside get_rna_base_centroid" );
		}
	} //if not RNA

	Vector centroid( 0.0 );
	Size numatoms = 0;

	//Consistency check:
	//if(rsd.type().atom_name(rsd.first_sidechain_atom()) !=" O2'") utility_exit_with_message( "rsd.type().atom_name(rsd.first_sidechain_atom()) !=\" O2'\" " );
	//if(rsd.atom_name( rsd.first_sidechain_atom() )!=" O2'") utility_exit_with_message("rsd.atom_name( rsd.first_sidechain_atom() )!=\" O2'\"");

	if ( rsd.RNA_info().o2prime_index() != rsd.first_sidechain_atom() ) {
		utility_exit_with_message( "rsd.RNA_info().o2prime_index() != rsd.first_sidechain_atom()" );
	}

	if ( verbose )  std::cout << "Base atoms" << std::endl;

	for ( Size i = rsd.first_sidechain_atom() + 1; i <= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.

		if ( verbose ) std::cout << "atom " << i  << " " <<  "name = " << rsd.type().atom_name( i ) << " type = " << rsd.atom_type( i ).name()  << " " << rsd.atom_type_index( i ) << " " << rsd.atomic_charge( i );

		if ( rsd.RNA_info().atom_is_virtual( i ) ) {
			if ( verbose ) std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}

		if ( rsd.atom_type( i ).is_repulsive() ) {
			if ( verbose ) std::cout << "  Repulsive: Ignore! " << std::endl;
			continue;
		}

		if ( verbose ) std::cout << std::endl;

		centroid += rsd.xyz( i );
		numatoms++;
	}

	if ( numatoms == 0 ) { //Centroid not well defined in this case...probably because rsd is a virtual residue...just return 0
		Vector dummy_centroid( 0.0 );
		return dummy_centroid;
	}

	centroid /= static_cast< Real > ( numatoms );

	return centroid;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Unify the version in StepWiseRNA_Util.cc and RNA_CentroidInfo.cc on June 25, 2011
numeric::xyzMatrix< core::Real >
get_rna_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ){

	using namespace chemical;

	//SML PHENIX conference
	if ( !rsd.is_RNA() ) {
		if ( basic::options::option[basic::options::OptionKeys::rna::erraser::rna_prot_erraser].value() ) {
			return numeric::xyzMatrix< core::Real > ::identity();
		} else { //if not option
			utility_exit_with_message( "non - RNA residue inside get_rna_base_coordinate_system, abort" );
		}
	} //if not RNA

	AA res_type = rsd.aa();
	if ( res_type == aa_unk || res_type == aa_unp ) res_type = rsd.type().na_analogue();

	Vector x, y, z;

	Vector WC_coord, H_coord;
	if ( res_type == aa_unk || res_type == aa_unp ) {
		// Just use the first two sidechain atoms for generality
		WC_coord = rsd.xyz( rsd.first_sidechain_atom() );
		H_coord = rsd.xyz( rsd.first_sidechain_atom() + 1 );
	} else {
		// Make an axis pointing from base centroid to Watson-Crick edge.
		std::string WC_atom;
		if ( res_type == na_rad || res_type == na_rgu ) WC_atom = "N1";
		if ( res_type == na_rcy || res_type == na_ura ) WC_atom = "N3";
		// Make a perpendicular axis pointing from centroid towards
		// Hoogstein edge (e.g., major groove in a double helix).
		std::string H_atom;
		if ( res_type == na_rad || res_type == na_rgu ) H_atom = "N7" ;
		if ( !rsd.has( H_atom ) && rsd.name3() == "7DA" ) H_atom = "C7";
		if ( res_type == na_rcy || res_type == na_ura ) H_atom = "C5";
		if ( res_type == na_ura && rsd.name3() == "PSU" ) H_atom = "N1"; // pseudoU is flipped

		WC_coord = rsd.xyz( WC_atom );
		H_coord = rsd.xyz( H_atom );
	}

	x = WC_coord - centroid;
	x.normalize();

	y = H_coord - centroid; //not orthonormal yet...
	z = cross( x, y );
	z.normalize(); // Should point roughly 5' to 3' if in a double helix.

	y = cross( z, x );
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
