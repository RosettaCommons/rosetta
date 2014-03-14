// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data_fixup.cc
/// @brief
/// @author Rhiju Das

// Unit headers
#include <core/io/pdb/file_data_fixup.hh>
#include <core/io/pdb/file_data.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// External headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <utility>

using namespace basic::options;
using namespace basic::options::OptionKeys;

// Tracer instance for this file
static basic::Tracer TR("core.io.pdb.file_data_fixup");

////////////////////////////////////////////////////////////////////////////////////////////
// useful utility scripts for fixing up residue and atom names.
//
// Originally implemented for nucleic acids, to permit backward-compatibility with
// old Rosetta atom & residue names.
//
// Some  useful functions below for fixing up H5'<-->H5'' ambiguities based on geometric
//  comparisons to ideal coordinates.
//
//
namespace core {
namespace io {
namespace pdb {

/// @details  Temporary hacky hack
/// Need better mechanism for this
/// for nucleic acids, slightly better mechanism is below (convert_nucleic_acid_residue_info_to_standard)
std::string
convert_res_name( std::string const & name )
{
	std::string name2 = name; //Copy of the input string for output.

	//Remove whitespace to make it easier to import metalloproteins:
	for(signed int i=0; i<(signed int)name2.length(); ++i) { //Needs to be signed!  Can't use core::Size!
		if(name2[i]==' ' || name2[i]=='\n') {
			name2.erase(i,1);
			--i;
		}
	}

	if ( name2 == "MSE" ) {
		TR << "Reading MSE as MET!" << std::endl;
		return "MET";
	}

	//If this is one of these metal ions and there is a 1 or 2 appended to the name (e.g. "CU2"), return just the first two characters as the name.
	//if(name2!="ZNx" && name2!="ZNX") { //Special zinc variants used in unit tests.  Grr.
	std::string const firsttwo = name2.substr(0,2);
	if (name2.length()>2 && (firsttwo == "CU" || firsttwo == "ZN" || firsttwo == "CA" || firsttwo == "CO" || firsttwo == "MG" || firsttwo == "MN" || firsttwo == "NA") ) {
		if(name2[2]=='1' || name2[2]=='2') {
			TR << "Reading " << name2 << " as " << firsttwo << "!" << std::endl;
			name2 = firsttwo;
		}
	}
	//}

	while(name2.size()<3) {
		name2 = std::string(" ") + name2; //Irritating -- name is expected to be exactly 3 characters.
	}

	return name2;
}

/// for nucleic acids, slightly better mechanism is below (convert_nucleic_acid_residue_info_to_standard)
std::string
convert_atom_name( std::string const & res_name, std::string atom_name )
{
	if( atom_name.size() != 4 ){
		std::string message= res_name+" has atom "+ atom_name+", with size!=4";
		utility_exit_with_message(message);
	};

	//atom_name = strip_whitespace( atom_name );
	if ( res_name == "MET" && atom_name == " S  " ) {
		return " SD ";
	} else if ( res_name == "MET" && atom_name == "SE  " ) {
		TR << "Reading Selenium SE from MSE as SD from MET" << std::endl;
		return " SD ";
	}
	return atom_name;
}

//////////////////////////////////////////////////////////////////////////////////////////
// @brief introduced in 2013 with changes of RNA and DNA atom types to match PDB -- rhiju.
//  Could also include cleanup/handling of chirality of hydrogens (e.g., H5' <--> H5'') in here,
//  but actually easier to do it later when we actually are instantiating ResidueType and have Rosetta's
//  ideal coordinates.
void
convert_nucleic_acid_residue_info_to_standard( 	utility::vector1< ResidueInformation > & all_rinfo, bool const force_RNA /* = false*/ ){

	// following is to show warnings or cap number.
	static Size nfix( 0 );
	static Size const max_fix( 2 );
	static bool const show_all_fixup( option[ in::show_all_fixes ]() );
	static bool showed_warning( false );

	for ( Size n = 1; n <= all_rinfo.size(); n++ ){

		ResidueInformation & rinfo  = all_rinfo[ n ];

		std::string const original_name = rinfo.resName;

		// first establish if this is DNA or RNA (or something else)
		if ( !force_RNA && is_potential_old_DNA( rinfo.resName ) && missing_O2prime( rinfo.atoms ) ) 	{
			rinfo.resName.replace( 1, 1, "D" ); // A --> dA
			if ( ++nfix <= max_fix || show_all_fixup ) TR << "Converting residue name " <<  original_name <<  " to " << rinfo.resName << std::endl;
			for ( Size n = 1 ; n <= rinfo.atoms.size(); n++ ) rinfo.atoms[ n ].resName = rinfo.resName;
		}

		if ( is_old_RNA( rinfo.resName ) ){
			rinfo.resName.replace( 1, 1, " " ); // rA --> A
			if ( ++nfix <= max_fix || show_all_fixup ) TR << "Converting residue name " <<  original_name <<  " to " << rinfo.resName << std::endl;
			for ( Size n = 1 ; n <= rinfo.atoms.size(); n++ ) rinfo.atoms[ n ].resName = rinfo.resName;
		}

		// Now apply atom mapping.
		if ( is_NA( rinfo.resName ) )	convert_nucleic_acid_atom_names_to_standard( rinfo, force_RNA );

	}

	if ( nfix > max_fix && !show_all_fixup && !showed_warning ){
		TR << "Number of nucleic acid residue fixups exceeds output limit. Rerun with -show_all_fixes to show everything." << std::endl;
		showed_warning = true;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////
bool
is_potential_old_DNA( std::string const & res_name ){
	return ( res_name == "  A" ||
					 res_name == "  C" ||
					 res_name == "  G" ||
					 res_name == "  T" );
}

//////////////////////////////////////////////////////////////////////////////////////////
bool
is_old_RNA( std::string const & res_name ){
	return ( res_name == " rA" ||
					 res_name == " rC" ||
					 res_name == " rG" ||
					 res_name == " rU" );
}

//////////////////////////////////////////////////////////////////////////////////////////
// for following, it might be smarter to define based on what is_NA() in residue type set (which
// is read from disk), rather than hard-wiring it.
bool
is_NA( std::string const & res_name ){
	return ( res_name == "  A" ||
					 res_name == "  C" ||
					 res_name == "  G" ||
					 res_name == "  U" ||
					 res_name == " DA" ||
					 res_name == " DC" ||
					 res_name == " DG" ||
					 res_name == " DT" );
}

//////////////////////////////////////////////////////////////////////////////////////////
bool
missing_O2prime( utility::vector1< AtomInformation > const & atoms ){

	for ( Size n = 1; n <= atoms.size(); n++ ){
		std::string const & name =  atoms[ n ].name;
		if ( name == " O2*" || name == " O2'" )  return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////
// @brief  This is a pretty good framework and could allow for other crazy nucleic acid atom name schemes.
void
convert_nucleic_acid_atom_names_to_standard( ResidueInformation & rinfo, bool const /* force_RNA=false */ )
{
	// following is to show warnings or cap number.
	static Size nfix( 0 );
	static Size const max_fix( 2 );
	static bool const show_all_fixup( option[ in::show_all_fixes ]() );
	static bool showed_warning( false );

	utility::vector1< AtomInformation > & all_atom_info = rinfo.atoms;

	for ( Size n = 1 ; n <= all_atom_info.size(); n++ ) {

		AtomInformation & atom_info = all_atom_info[n];

		std::string const original_name = atom_info.name;

		convert_nucleic_acid_atom_name_to_standard( atom_info );

		std::string const new_atom_name = atom_info.name;
		if ( original_name != new_atom_name ){
			if ( ++nfix <= max_fix || show_all_fixup ) TR << "Converting atom name    " << original_name << " to " << new_atom_name << std::endl;

			rinfo.xyz[ new_atom_name ]   = rinfo.xyz[ original_name ];
			rinfo.xyz.erase( original_name );

			rinfo.temps[ new_atom_name ] = rinfo.temps[ original_name ];
			rinfo.temps.erase( original_name );
		}

	}

	if ( nfix > max_fix && !show_all_fixup && !showed_warning ){
		TR << "Number of atom_name fixups exceeds output limit. Rerun with -show_all_fixes to show everything." << std::endl;
		showed_warning = true;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////
// following could also be a class...
void
convert_nucleic_acid_atom_name_to_standard( AtomInformation & atom_info ){

	static std::map < std::string, std::string > atom_name_map;
	static bool init( false );


	//  stars (*)  are changed to primes (').
	if ( atom_info.name[3] == '*' ) atom_info.name = atom_info.name.substr(0,3) + "\'";

	// just initialize map once.
	if ( ! init ){

		atom_name_map[ " O1P" ] = " OP2";
		atom_name_map[ " O2P" ] = " OP1";
		atom_name_map[ "1H2'" ] = " H2'";
		atom_name_map[ "1H5'" ] = " H5'";
		atom_name_map[ "2H5'" ] = "H5''";
		atom_name_map[ "2HO'" ] = "HO2'";

		atom_name_map[ "1H4 " ] = " H41";
		atom_name_map[ "2H4 " ] = " H42";
		atom_name_map[ "1H6 " ] = " H61";
		atom_name_map[ "2H6 " ] = " H62";

		atom_name_map[ "2H2 " ] = " H21";
		atom_name_map[ "1H2 " ] = " H22";

		// note this will end up with the wrong chirality
		// for DNA, but we'll fix that later.
		atom_name_map[ "2H2'" ] = "H2''";

		// thymidine.
		atom_name_map[ " C5M" ] = " C7 ";
		atom_name_map[ "1H5M" ] = " H71";
		atom_name_map[ "2H5M" ] = " H72";
		atom_name_map[ "3H5M" ] = " H73";

		init = true;
	}

	if ( atom_name_map.find( atom_info.name ) != atom_name_map.end() ) atom_info.name = atom_name_map[ atom_info.name ];

}



//////////////////////////////////////////////////////////////////////////////////////////
// @brief due to differences in different crystallography/NMR/modeling packages, labeling of sister atoms
//  (like OP1 <--> OP2, or H41 <--> H42) in PDBs is totally wacky. This is an attempt to regularize...
//  and it can actually make a difference since sometimes partial charges on sister hydrogens
//  can be different. Right now only set up for nucleic acids, but could probably generalize.
void
check_and_correct_sister_atoms( core::conformation::ResidueOP & rsd ){

	// in all nucleic acids
	check_and_correct_sister_atom_based_on_chirality( rsd, " OP1", " OP2", " P  ", " O5'" );
	check_and_correct_sister_atom_based_on_chirality( rsd, " H5'", "H5''", " C5'", " C4'" );
	// in DNA
	check_and_correct_sister_atom_based_on_chirality( rsd, " H2'", "H2''", " C2'", " C3'" );

	// in adenosine
	check_and_correct_sister_atom_based_on_outgroup( rsd, " H61", " H62", " N1 " );
	// in guanosine
	check_and_correct_sister_atom_based_on_outgroup( rsd, " H21", " H22", " N1 " );
	// in cytidine
	check_and_correct_sister_atom_based_on_outgroup( rsd, " H41", " H42", " N3 " );

}


//////////////////////////////////////////////////////////////////////////////////////////
// sisters sprout off the same parent, and outer_ref is something else bonded to the parent. a cousin, i guess.
void
check_and_correct_sister_atom_based_on_chirality( core::conformation::ResidueOP & rsd,
																									std::string const sister1_name,
																									std::string const sister2_name,
																									std::string const parent_name,
																									std::string const cousin_name ){

	if ( !rsd->has( sister1_name ) ) return;
	if ( !rsd->has( sister2_name ) ) return;
	if ( !rsd->has( parent_name ) ) return;
	if ( !rsd->has( cousin_name ) ) return;

 	Vector const current_xyz_sister1        = rsd->xyz( sister1_name );
	Vector const current_xyz_sister2        = rsd->xyz( sister2_name );
	Vector const current_xyz_parent         = rsd->xyz( parent_name );
	Vector const current_xyz_cousin      = rsd->xyz( cousin_name );
	int current_sign = get_chirality_sign( current_xyz_sister1, current_xyz_sister2, current_xyz_parent, current_xyz_cousin );

	core::chemical::ResidueType const & rsd_type = rsd->type();
 	Vector const ideal_xyz_sister1        = rsd_type.atom( sister1_name ).ideal_xyz();
	Vector const ideal_xyz_sister2        = rsd_type.atom( sister2_name ).ideal_xyz();
	Vector const ideal_xyz_parent         = rsd_type.atom( parent_name ).ideal_xyz();
	Vector const ideal_xyz_cousin      = rsd_type.atom( cousin_name ).ideal_xyz();
	int ideal_sign = get_chirality_sign( ideal_xyz_sister1, ideal_xyz_sister2, ideal_xyz_parent, ideal_xyz_cousin );

	if ( current_sign != ideal_sign )	flip_atom_xyz( rsd, sister1_name, sister2_name );

}

//////////////////////////////////////////////////////////////////////////////////
void
check_and_correct_sister_atom_based_on_outgroup( core::conformation::ResidueOP & rsd,
																								 std::string const sister1_name,
																								 std::string const sister2_name,
																								 std::string const outgroup_name ){

	if ( !rsd->has( sister1_name ) ) return;
	if ( !rsd->has( sister2_name ) ) return;
	if ( !rsd->has( outgroup_name ) ) return;

 	Vector const current_xyz_sister1        = rsd->xyz( sister1_name );
	Vector const current_xyz_sister2        = rsd->xyz( sister2_name );
	Vector const current_xyz_outgroup       = rsd->xyz( outgroup_name );

	int current_closest_sister = get_closest_sister( current_xyz_sister1, current_xyz_sister2, current_xyz_outgroup );

	core::chemical::ResidueType const & rsd_type = rsd->type();
 	Vector const ideal_xyz_sister1        = rsd_type.atom( sister1_name ).ideal_xyz();
	Vector const ideal_xyz_sister2        = rsd_type.atom( sister2_name ).ideal_xyz();
	Vector const ideal_xyz_outgroup       = rsd_type.atom( outgroup_name ).ideal_xyz();
	int ideal_closest_sister = get_closest_sister( ideal_xyz_sister1, ideal_xyz_sister2, ideal_xyz_outgroup );

	if ( current_closest_sister != ideal_closest_sister )	flip_atom_xyz( rsd, sister1_name, sister2_name );

}

//////////////////////////////////////////////////////////////////////////////////
void
flip_atom_xyz( core::conformation::ResidueOP & rsd,
							 std::string const & sister1_name,
							 std::string const & sister2_name ) {
	// following is to show warnings or cap number.
	static Size nfix( 0 );
	static Size const max_fix( 2 );
	static bool const show_all_fixup( option[ in::show_all_fixes ]() );
	static bool showed_warning( false );

	if ( ++nfix <= max_fix || show_all_fixup )	{
		TR << "Flipping atom xyz for " << sister1_name << " and " << sister2_name << " for residue " << rsd->name3() <<  std::endl;
	}
	Vector const temp_xyz = rsd->xyz( sister1_name );
	rsd->set_xyz( sister1_name, rsd->xyz( sister2_name ) );
	rsd->set_xyz( sister2_name, temp_xyz );

	if ( nfix > max_fix && !show_all_fixup && !showed_warning ){
		TR << "Number of flip-atom fixups exceeds output limit. Rerun with -show_all_fixes to show everything." << std::endl;
		showed_warning = true;
	}
}

//////////////////////////////////////////////////////////////////////////////////
int
sgn( Real const & x ){
	return ( x > 0 ) - ( x < 0 );
}

//////////////////////////////////////////////////////////////////////////////////
int
get_chirality_sign(  Vector const & xyz_sister1,
										 Vector const & xyz_sister2,
										 Vector const & xyz_parent,
										 Vector const & xyz_cousin ) {
	int const sign = sgn( dot( xyz_cousin - xyz_parent, cross( xyz_sister1 - xyz_parent, xyz_sister2 - xyz_parent ) ) );
	if ( sign == 0 ) utility_exit_with_message( "unexpected sign error when checking chirality" );
	return sign;
}

//////////////////////////////////////////////////////////////////////////////////
// returns 1 or 2 based on which sister is closest to outgroup.
int
get_closest_sister(  Vector const & xyz_sister1,
									 Vector const & xyz_sister2,
									 Vector const & xyz_outgroup ) {
	return ( xyz_sister1.distance( xyz_outgroup ) < xyz_sister2.distance( xyz_outgroup) ) ? 1 : 2;
}

} // namespace pdb
} // namespace io
} // namespace core
