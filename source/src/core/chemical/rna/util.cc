// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rna/RNA_Util.cc
/// @author Rhiju Das

// Unit headers
#include <core/chemical/rna/util.hh>
#include <core/types.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>

// Project headers
#include <numeric/constants.hh>

#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

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
	if ( c == 'a' ) return 1;
	if ( c == 'c' ) return 2;
	if ( c == 'g' ) return 3;
	if ( c == 'u' ) return 4;
	return 0;
}

/////////////////////////////////////////////////////////////////////////////
//This may be used elsewhere -- set up a util.hh?
char get_edge_from_num( Size const num ) {
  if ( num == WATSON_CRICK ) return 'W';
  if ( num == HOOGSTEEN )    return 'H';
  if ( num == SUGAR )        return 'S';
  if ( num == O2PRIME )       return '2';
  if ( num == PHOSPHATE )    return 'P';
  return 'X';
}

/////////////////////////////////////////////////////////////////////////////
std::string  //Parin March 7, 2011
get_full_edge_from_num( Size const num ) {
  if ( num == WATSON_CRICK ) return "WC";
  if ( num == HOOGSTEEN )    return "HOOG";
  if ( num == SUGAR )        return "SUGAR";
  if ( num == O2PRIME )       return "O2PRIME";
  if ( num == PHOSPHATE )    return "PHOS";

	std::cout << "Invalid edge num = " << num << std::endl;
	utility_exit_with_message( "Invalid edge num!" );
  return "ERROR";

}

/////////////////////////////////////////////////////////////////////////////
//This may be used elsewhere -- set up a util.hh?
char get_orientation_from_num( Size const num ) {
	if ( num == 1 ) return 'A';
	if ( num == 2 ) return 'P';
	return 'X';
}

/////////////////////////////////////////////////////////////////////////////

std::string //Parin March 7, 2011
get_full_orientation_from_num( Size const num ) {
  if ( num == 0 ) return "BLAH";
  if ( num == 1 ) return "ANTI";
  if ( num == 2 ) return "PARA";

	std::cout << "Invalid orientation num = " << num << std::endl;
	utility_exit_with_message( "Invalid orientation num!" );
  return "ERROR";
}

std::string //Parin April 19, 2011
get_full_LW_orientation_from_num( Size const num ){
  if ( num == 0 ) return "BLAH ";
  if ( num == 1 ) return "CIS  ";
  if ( num == 2 ) return "TRANS";

	std::cout << "Invalid orientation num = " << num << std::endl;
	utility_exit_with_message( "Invalid orientation num!" );
  return "ERROR";
}

///////////////////////////////////////////////////////////////////////////////
std::string const	first_base_atom( conformation::Residue const & rsd ) {
	//	if (rsd.name1() == 'a' || rsd.name1() == 'g' ) 	return " N9 ";
	//	return " N1 ";
	return rsd.atom_name( first_base_atom_index( rsd ) );
}

///////////////////////////////////////////////////////////////////////////////
bool	is_purine( conformation::Residue const & rsd ) {
	if ( rsd.name1() == 'a' || rsd.name1() == 'g' ) 	return true;
	return false;
}


///////////////////////////////////////////////////////////////////////////////
Size first_base_atom_index( conformation::Residue const & rsd ) {
	// HEY MAKE THIS MORE GENERAL? Maybe look at chi1?
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 1 /*chi # 1 must be nucleic acid "chi"*/ );
	return atom_indices[ 3 ]; /* C2' ... C1' ... first base atom ...  chi1 torsion atom*/
}

///////////////////////////////////////////////////////////////////////////////
std::string const	chi1_torsion_atom( conformation::Residue const & rsd ) {
	//	if (rsd.name1() == 'a' || rsd.name1() == 'g' ) 	return " N9 ";
	//	return " N1 ";
	return rsd.atom_name( chi1_torsion_atom_index( rsd ) );
}


///////////////////////////////////////////////////////////////////////////////
Size chi1_torsion_atom_index( conformation::Residue const & rsd ) {
	// HEY MAKE THIS MORE GENERAL? Maybe look at chi1?
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 1 /*chi # 1 must be nucleic acid "chi"*/ );
	return atom_indices[ 4 ]; /* C2' ... C1' ... first base atom ... chi1 torsion atom*/
}


///////////////////////////////////////////////////////////////////////////////
// consider moving this to chemical/util.cc.
std::string const	default_jump_atom( conformation::Residue const & rsd ) {
	if ( rsd.is_RNA() ){
		if ( !rsd.is_coarse() ){
			return chi1_torsion_atom( rsd );
		} else {
			return " Y  ";
		}
	}
	if ( rsd.is_protein() ){
		return " CA "; // note that this does not match 'traditional' choice in FoldTree.cc
	}
	if ( rsd.name3() == " MG" )		return "MG  ";
	if ( rsd.name3() == " ZN" )		return "ZN  ";
	if ( rsd.name3() == "XXX" )	return " Y  ";

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

///////////////////////////////////////////////////////////////////////////////
void
get_watson_crick_base_pair_atoms(
	 chemical::AA const & aa1,
	 chemical::AA const & aa2,
	 std::string & atom1,
	 std::string & atom2 ) {

	using namespace core::chemical;

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
	return;
}

/////////////////////////////////////////////////////////////////////
void
get_watson_crick_base_pair_atoms(
	 chemical::AA const & aa1,
	 chemical::AA const & aa2,
	 utility::vector1< std::string > & atom_ids1,
	 utility::vector1< std::string > & atom_ids2	 )
{

	using namespace chemical;

	atom_ids1.clear();
	atom_ids2.clear();

	if ( aa1 == na_rad && aa2 == na_ura ) {
		atom_ids1.push_back( " N1 " );		atom_ids2.push_back( " H3 " );
		atom_ids1.push_back( " H61" );		atom_ids2.push_back( " O4 " );
		return;
	} else if ( aa1 == na_rgu && aa2 == na_rcy ) {
		atom_ids1.push_back( " H1 " );		atom_ids2.push_back( " N3 " );
		atom_ids1.push_back( " H21" );		atom_ids2.push_back( " O2 " );
		atom_ids1.push_back( " O6 " );		atom_ids2.push_back( " H41" );
		return;
	} else if ( aa1 == na_rgu && aa2 == na_ura ) {
		atom_ids1.push_back( " O6 " );		atom_ids2.push_back( " H3 " );
		atom_ids1.push_back( " H1 " );		atom_ids2.push_back( " O2 " );
		return;
	} else	if ( aa2 == na_rad && aa1 == na_ura ) {
		atom_ids2.push_back( " N1 " );		atom_ids1.push_back( " H3 " );
		atom_ids2.push_back( " H61" );		atom_ids1.push_back( " O4 " );
		return;
	} else if ( aa2 == na_rgu && aa1 == na_rcy ) {
		atom_ids2.push_back( " H1 " );		atom_ids1.push_back( " N3 " );
		atom_ids2.push_back( " H21" );		atom_ids1.push_back( " O2 " );
		atom_ids2.push_back( " O6 " );		atom_ids1.push_back( " H41" );
		return;
	} else if ( aa2 == na_rgu && aa1 == na_ura ) {
		atom_ids2.push_back( " O6 " );		atom_ids1.push_back( " H3 " );
		atom_ids2.push_back( " H1 " );		atom_ids1.push_back( " O2 " );
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
		if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
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

	if ( rsd.RNA_type().o2prime_index() != rsd.first_sidechain_atom() ){
		utility_exit_with_message( "rsd.RNA_info().o2prime_index() != rsd.first_sidechain_atom()" );
	}

	if ( verbose )  std::cout << "Base atoms" << std::endl;

	for ( Size i = rsd.first_sidechain_atom() + 1; i <= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.

		if ( verbose ) std::cout << "atom " << i  << " " << 	"name = " << rsd.type().atom_name( i ) << " type = " << rsd.atom_type( i ).name()  << " " << rsd.atom_type_index( i ) << " " << rsd.atomic_charge( i );

		if ( rsd.RNA_type().atom_is_virtual( i ) ){
			if ( verbose ) std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}

		if ( verbose ) std::cout << std::endl;

    centroid += rsd.xyz( i );
    numatoms++;
  }

	if ( numatoms == 0 ){ //Centroid not well defined in this case...probably because rsd is a virtual residue...just return 0
		Vector dummy_centroid( 0.0 );
		return dummy_centroid;
	}

  centroid /= static_cast< Real > ( numatoms );

  return centroid;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Unify the version in StepWiseRNA_Utill.cc and RNA_CentroidInfo.cc on June 25, 2011
numeric::xyzMatrix< core::Real >
get_rna_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ){

	using namespace chemical;

  //SML PHENIX conference
  if ( !rsd.is_RNA() ) {
    if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
      return numeric::xyzMatrix< core::Real > ::identity();
    } else { //if not option
      utility_exit_with_message( "non - RNA residue inside get_rna_base_coordinate_system, abort" );
    }
  } //if not RNA

 	Size res_type = rsd.aa();

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
		if ( res_type == na_rcy || res_type == na_ura ) H_atom = "C5";
		WC_coord = rsd.xyz( WC_atom );
		H_coord = rsd.xyz( H_atom );
	}

	x = WC_coord - centroid;
	x.normalize();

	y = H_coord - centroid; //not orthonormal yet...
	z = cross( x, y );
	z.normalize(); // Should poSize roughly 5' to 3' if in a double helix.

	y = cross( z, x );
	y.normalize(); //not necessary but doesn't hurt.

 	numeric::xyzMatrix< core::Real > M = numeric::xyzMatrix< core::Real > ::cols( x, y, z );
	return M;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Check whether one of the atom beyond to a base and another to a phosphate. Doesn't have to be on the same nucleotide.
bool
Is_base_phosphate_atom_pair( conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const atomno_1, Size const atomno_2 ){

	bool Is_base_phosphate_atom_pair = false;

	if ( ( rsd_1.RNA_type().atom_is_phosphate( atomno_1 ) && ( rsd_2.RNA_type().is_RNA_base_atom( atomno_2 ) ) ) ) Is_base_phosphate_atom_pair = true;
	if ( ( rsd_2.RNA_type().atom_is_phosphate( atomno_2 ) && ( rsd_1.RNA_type().is_RNA_base_atom( atomno_1 ) ) ) ) Is_base_phosphate_atom_pair = true;

	if ( Is_base_phosphate_atom_pair ){ //This Assume that rsd_1 and rsd_2 are the same!!!
		if ( rsd_1.seqpos() == rsd_2.seqpos() && ( rsd_1.path_distance( atomno_1, atomno_2 ) < 4 ) ){ //consistency check!
			utility_exit_with_message( "Is_base_phosphate_atom_pair but rsd.path_distance( " + ObjexxFCL::string_of( atomno_1 ) + ", " + ObjexxFCL::string_of( atomno_2 ) + " ) < 4" );
		}
	}

	return Is_base_phosphate_atom_pair;

}


} //ns rna
} //ns chemical
} //ns core
