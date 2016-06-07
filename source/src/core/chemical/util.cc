// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/util.cc
/// @brief Utilities for modifying and utilizing ResidueTypes and other core::chemical classes.

// Unit headers
#include <core/chemical/util.hh>

// Package Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/Patch.hh>

// Project Headers
#include <core/types.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/string.functions.hh>


namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.util" );

// Return a constant access pointer to the ResidueTypeSet specified by the command-line options.
core::chemical::ResidueTypeSetCAP
rsd_set_from_cmd_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const type_set_name( option[ in::file::residue_type_set ]() );
	ResidueTypeSetCAP set = ChemicalManager::get_instance()->residue_type_set( type_set_name );

	return set;
}


// Add additional parameter files not present in <atom-set-name>/extras.txt.
/// @details Called by ChemicalManager at time of AtomTypeSet creation.
void
add_atom_type_set_parameters_from_command_line(
	std::string const & atom_type_set_tag,
	AtomTypeSet & atom_type_set )
{
	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::chemical;

	if ( ! ( option[ add_atom_type_set_parameters ].user() ) ) {
		return;  // do nothing if flag not present
	}
	utility::vector1< std::string > paramstring( option[ add_atom_type_set_parameters ]() );
	if ( paramstring.size() % 2 != 0 ) {
		utility_exit_with_message( "bad format for -add_atom_type_set_parameters; "
			"should be: -add_atom_type_set_parameters <tag1> <filename1> <tag2> <filename2> ...");
	}
	Size const nparams( paramstring.size() / 2 );
	for ( core::uint i = 0; i < nparams; ++i ) {
		string const tag( paramstring[ 2 * i + 1 ] ), filename( paramstring[ 2 * i + 2 ] );
		TR.Trace << "add_atom_type_set_parameters_from_command_line: desired_tag= " << atom_type_set_tag <<
			" cmdline-tag= " << tag << " filename= " << filename << endl;
		if ( tag == atom_type_set_tag ) {
			if ( ! utility::file::file_exists( filename ) ) {
				utility_exit_with_message( "unable to locate/open file: " + filename );
			}
			TR.Trace << "add_atom_type_set_parameters_from_command_line: tag= " << tag << " filename= " << filename <<
				endl;
			atom_type_set.add_parameters_from_file( filename );
		}
	}
}


// Modify atom_type properties from the command line.
/// @details  Called by ChemicalManager at time of AtomTypeSet creation.
void
modify_atom_properties_from_command_line(
	std::string const & atom_type_set_tag,
	AtomTypeSet & atom_type_set )
{
	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::chemical;
	using namespace ObjexxFCL;

	if ( option[ set_atom_properties ].user() ) {
		utility::vector1< std::string > const & mods( option[ set_atom_properties ] );

		string const errmsg( "-set_atom_properties format should be:: -set_atom_properties <set1>:<atom1>:<param1>:"
			"<setting1> <set2>:<atom2>:<param2>:<setting2> ...; for example: '-chemical:set_atom_properties "
			"fa_standard:OOC:LK_DGFREE:-5 fa_standard:ONH2:LJ_RADIUS:0.5' ");

		for ( core::uint i = 1; i <= mods.size(); ++i ) {
			// mod should look like (for example):  "fa_standard:OOC:LK_RADIUS:4.5"
			string const & mod( mods[ i ] );

			core::uint const pos1( mod.find( ":" ) );
			if ( pos1 == string::npos ) utility_exit_with_message( errmsg );
			string const atomset_tag( mod.substr( 0, pos1 ) );
			if ( atomset_tag != atom_type_set_tag ) continue;

			core::uint const pos2( mod.substr( pos1 + 1 ).find( ":" ) );
			if ( pos2 == string::npos ) utility_exit_with_message( errmsg );
			string const atom_name( mod.substr( pos1 + 1, pos2 ) );
			if ( ! atom_type_set.has_atom( atom_name ) ) {
				utility_exit_with_message( errmsg + ". Nonexistent atomname: " + atom_name );
			}

			core::uint const pos3( mod.substr( pos1 + 1 + pos2 + 1 ).find( ":" ) );
			if ( pos3 == string::npos ) utility_exit_with_message( errmsg );
			string const param( mod.substr( pos1 + 1 + pos2 + 1, pos3 ) );

			string const stringsetting( mod.substr( pos1 + 1 + pos2 + 1 + pos3 + 1 ) );
			if ( ! is_double( stringsetting ) ) utility_exit_with_message( errmsg );
			Real const setting( double_of( stringsetting ) );

			TR.Trace << "modify_atom_properties_from_command_line: setting " << atomset_tag << ' ' << atom_name <<
				' ' << param << ' ' << setting << endl;

			Size const atom_index( atom_type_set.atom_type_index( atom_name ) );

			// I would like to uncomment the following if-check, but right now there is an extra parameter file
			// that defines a parameter with the name LK_DGFREE (memb_fa_params.txt). That's kind of confusing...
			//
			// if ( atom_type_set.has_extra_parameter( param ) ) {
			//  Size const param_index( atom_type_set.extra_parameter_index( param ) );
			//  atom_type_set[ atom_index ].set_extra_parameter( param_index, setting );
			// } else {
			atom_type_set[ atom_index ].set_parameter( param, setting );
		}
	}
}


// Return a string representing the internal coordinates tree of this ResidueType.
/// @note  Mainly intended for debugging purposes.
std::string
formatted_icoord_tree( core::chemical::ResidueType const & restype )
{
	using namespace std;
	using namespace utility;

	// General scheme: pick an atom, then move up the icoord tree until you hit the root.
	// Then proceed down the tree depth first, outputting atom names as you go, keeping a stack of atoms to go back to.
	// The complicated bit is keeping track of where the different attachment points are, such that you can accurately
	// place the parenthesis for the tree.
	string output( "Icoord tree: " );

	vector1< bool > placed( restype.natoms(), false );
	vector1< AtomICoor > icoords;
	vector1< core::Size > deferred;
	for ( core::uint ii( 1 ); ii <= restype.natoms(); ++ii ) {
		icoords.push_back( restype.icoor( ii ) );
	}

	core::uint curres( 1 );  // The first atom is arbitrary (although it probably is the root).
	// The root atom has itself listed as the stub_atom1.
	while ( icoords[ curres ].stub_atom1().atomno() != curres ) {
		curres = icoords[ curres ].stub_atom1().atomno();
	}

	do {
		output += restype.atom_name( curres );
		placed[ curres ] = true;

		vector1< Size > possibles;
		for ( core::uint ii( 1 ); ii <= icoords.size(); ++ii ) {
			if ( icoords[ ii ].stub_atom1().atomno() == curres && ! placed[ ii ] ) {
				possibles.push_back( ii );
			}
		}

		if ( possibles.size() == 1 ) {
			curres = possibles[ 1 ];
		} else if ( possibles.size() ) {
			output += " (";
			curres = possibles.back();
			possibles.pop_back();
			if ( deferred.size() ) {
				deferred.push_back( 0 ); // Sentinel, but not needed for first
			}
			deferred.insert( deferred.end() , possibles.begin(), possibles.end() );
		} else if ( deferred.size() ) {
			output += ") ";
			while ( deferred.size() && deferred.back() == 0 ) {
				output += ") ";
				deferred.pop_back();
			}
			if ( ! deferred.size() ) {
				curres = 0;
				break;
			}
			output += " (";
			curres = deferred.back();
			deferred.pop_back();
		} else {
			output += ") ";
			curres = 0;
			break;
		}
	} while( curres != 0 );

	return output;
}


// Utility to examine chi output.
void
print_chis( std::ostream & out, ResidueType const & res )
{
	using namespace std;
	using namespace core::chemical;

	out << "Residue: " << res.name() << endl;
	out << formatted_icoord_tree(res) << endl;
	for ( core::uint ii( 1 ); ii <= res.nchi(); ++ii ) {
		AtomIndices const & indexes( res.chi_atoms( ii ) );
		out << "Chi " << ii << ": " << res.atom_name( indexes[ 1 ] ) << " " << res.atom_name( indexes[ 2 ] ) << " "
			<< res.atom_name( indexes[ 3 ] ) << " " << res.atom_name( indexes[ 4 ] ) << " ";
		if ( res.is_proton_chi( ii ) ) { out << " PROTON"; }
		out << endl;
	}
}


// Replaces the deprecated "_p:" linker connecting ResidueType base names with their patch names with ":".
/// @note This is here for backwards compatibility.
std::string
fixup_patches( std::string string_in )
{
	std::string string_out = string_in;
	string_out = utility::replace_in( string_out, "_p:", PATCH_LINKER );
	string_out = utility::replace_in( string_out, "Virtual_RNA_Residue_Upper", "Virtual_Phosphate" );
	string_out = utility::replace_in( string_out, "Virtual_Phosphate"+PATCH_LINKER+"Virtual_Phosphate", "Virtual_Phosphate" );
	return string_out;
}


// Are these two residues patched in exactly the same way, ignoring any VariantTypes in the list of exceptions?
/// @author  Labonte <JWLabonte@jhu.edu>
bool
variants_match_with_exceptions(
	ResidueType const & res1,
	ResidueType const & res2,
	utility::vector1< VariantType > list_of_variants_to_ignore )
{
	using namespace std;
	using namespace utility;

	for ( VariantType variant = FIRST_VARIANT; variant <= N_VARIANTS; ++variant ) {
		if ( list_of_variants_to_ignore.has_value( variant ) ) continue;
		if ( res1.properties().is_variant_type( variant ) != res2.properties().is_variant_type( variant ) ) {
			return false;
		}
	}

	// Now check for "custom" VariantTypes, which are stored as strings.
	if ( res1.properties().has_custom_variant_types() || res2.properties().has_custom_variant_types() ) {
		vector1 < string > const & list1( res1.properties().get_list_of_custom_variants() ),
			list2( res2.properties().get_list_of_custom_variants() );
		Size const n_variants_in_list1( list1.size() ),
			n_variants_in_list2( list2.size() );
		if ( n_variants_in_list1 != n_variants_in_list2 ) {
			return false;
		} else {
			for ( core::uint i = 1; i <= n_variants_in_list1; ++i ) {
				if ( ! list2.has_value( list1[ i ] ) ) {
					return false;
				}
			}
		}
	}

	return true;
}

/// @brief   check if user has set -pH_mode.
/// @details used to determine exceptions to PROTONATION/DEPROTONAT in variants_match.
utility::vector1< VariantType >
pH_mode_exceptions() {
	utility::vector1< VariantType > exceptions;
	using namespace basic::options;
	if ( option[ OptionKeys::pH::pH_mode ]() ) {
		exceptions.push_back( PROTONATED );
		exceptions.push_back( DEPROTONATED );
	}
	return exceptions;
}

// Are these two residues patched in exactly the same way?
/// @details If pH mode is being used, this function returns true, even if the two residues compared have different
/// protonation states.
bool
variants_match( ResidueType const & res1, ResidueType const & res2 )
{
	return variants_match_with_exceptions( res1, res2, pH_mode_exceptions() );
}


// Similar to variants_match(), but allows different adduct-modified states.
bool
nonadduct_variants_match( ResidueType const & res1, ResidueType const & res2 )
{
	return variants_match_with_exceptions( res1, res2, utility::vector1< VariantType >( 1, ADDUCT_VARIANT ) );
}



///////////////////////////////////////////////////////////////////////////////
///
/// @brief look for best match to atom_names
/// @details taken out of build_pose_as_is1
///  rsd_type should have all the atoms present in xyz
///  try to minimize atoms missing from xyz
ResidueTypeCOP
find_best_match( ResidueTypeCOPs const & rsd_type_list,
	utility::vector1< std::string > const & atom_names,
	bool const ignore_atom_named_H /* = false */ )
{

	using namespace core::chemical;
	Size best_index(0), best_rsd_missing( 99999 ), best_xyz_missing( 99999 );
	for ( Size j=1; j <= rsd_type_list.size(); j++ ) {

		ResidueType const & rsd_type( *(rsd_type_list[j]) );

		Size rsd_missing( 0 ), xyz_missing( 0 );

		for ( Size k=1; k<= rsd_type.natoms(); ++k ) {
			bool found_match( false );
			for ( Size m = 1; m <= atom_names.size(); ++m ) {
				if ( ObjexxFCL::stripped_whitespace( atom_names[m] ) == ObjexxFCL::stripped_whitespace( rsd_type.atom_name(k) ) ) {
					found_match = true; break;
				}
			}
			if ( !found_match ) ++xyz_missing;
		}

		for ( Size n = 1; n <= atom_names.size(); n++ ) {
			std::string const & atom_name = atom_names[ n ];
			if ( !rsd_type.has( ObjexxFCL::stripped_whitespace( atom_name ) ) &&
					!( atom_name == " H  " && ignore_atom_named_H ) ) { // don't worry about missing BB H if Nterm
				++rsd_missing;
			}
		}

		if ( ( rsd_missing < best_rsd_missing ) ||
				( rsd_missing == best_rsd_missing && xyz_missing < best_xyz_missing ) ) {
			best_rsd_missing = rsd_missing;
			best_xyz_missing = xyz_missing;
			best_index = j;
		}
	} // j=1,rsd_type_list.size()

	return  rsd_type_list[ best_index ];
}

//////////////////////////////////////
// rhiju/fang -- Use larger LJ_WDEPTH for protons to avoid clashes in RNA
void
enlarge_h_lj_wdepth( utility::vector1< Real > & lj_wdepth, AtomTypeSet const & atom_type_set ) {
	Real const enlarged_lj_wdepth = 0.15;
	// why not just look at element type? Anyway...
	utility::vector1< std::string > const H_names = utility::tools::make_vector1( "Hpol", "Hapo", "Haro", "HNbb", "HOH" );
	runtime_assert( lj_wdepth.size() == atom_type_set.n_atomtypes() );
	for ( Size i = 1; i <= H_names.size(); ++i ) {
		Size const index = atom_type_set.atom_type_index( H_names[i] );
		lj_wdepth[ index ] = enlarged_lj_wdepth;
	}
}

//////////////////////////////////////
// rhiju/fang -- Use larger LJ_WDEPTH for protons to avoid clashes in RNA
// this is a wrapper around another function to avoid copying code, including
// Fang's choice of the enlarged_lj_wdepth parameter.
void
enlarge_h_lj_wdepth( AtomTypeSet & atom_type_set ) {
	utility::vector1< Real > lj_wdepth;
	for ( Size n = 1; n <= atom_type_set.n_atomtypes(); n++ ) lj_wdepth.push_back( atom_type_set[n].lj_wdepth() );
	enlarge_h_lj_wdepth( lj_wdepth, atom_type_set );
	for ( Size n = 1; n <= atom_type_set.n_atomtypes(); n++ ) atom_type_set[n].set_parameter( "LJ_WDEPTH", lj_wdepth[n] );
}


//////////////////////////////////////
// In RNA modeling, getting unusual H-bonds to O4', O5', and O3' ...
// those are ethers and probably not great hydrogen bond acceptors (except maybe in
// catalytic transition states)
void
turn_off_hbonds_to_ether_oxygens( AtomTypeSet & atom_type_set ) {
	utility::vector1< Real > lj_wdepth;
	utility::vector1< std::string > const Oet_names = utility::tools::make_vector1( "Oet2", "Oet3" );
	std::string property; // have to set through a string since set_property() accepts property by reference.
	for ( Size i = 1; i <= Oet_names.size(); ++i ) {
		Size const index = atom_type_set.atom_type_index( Oet_names[i] );
		property = "ACCEPTOR";  atom_type_set[ index ].set_property( property, false );
		property = "DONOR"   ;  atom_type_set[ index ].set_property( property, false );
	}
}

void
detect_ld_chirality_from_polymer_residue(
	std::map< std::string, Vector > const & xyz,
	std::string const & name3,
	bool & is_d_aa,
	bool & is_l_aa
) {
	// We need a list of name3s with a given property like RNA.
	// So basically you could ask if a string has a given property
	// this way we don't have to have these disgusting things
	// We would populate this from the database.

	// Also, we need this function not to call on everything including e.g. water,
	// ions...

	is_d_aa = false;
	is_l_aa = false;

	// Exclude known achiral.
	if ( name3 == "GLY" || name3 == "C15" || name3 == "C16" || name3 == "MAL" ||
			name3 == "A98" || name3 == "B02" || name3 == "B06" ) {
		return;
	}

	// return false,false for peptoids and pna
	if ( xyz.find( " NG " ) != xyz.end() && ( name3 == "UPN" || name3 == "APN" || name3 == "TPN" || name3 == "GPN" || name3 == "CPN" ) ) {
		return;
	}

	// If termini are missing, then properly speaking chirality can't be inferred at all.
	// Let this be an unrecognized residue (as before).

	// Note that because there is that one phosphonate-Cterm residue
	// we have to check for EITHER Pbb or C being there.

	// AMW: no, we have to assume L or we break chainbreak stuff?

	// OK, we need ONE or more but not three protein bb atoms.
	// This will protect against carbohydrates and nucleic acids.
	Size bb_atoms_found = 0;
	if ( xyz.find( " N  " ) != xyz.end() ) ++bb_atoms_found;
	if ( xyz.find( " CA " ) != xyz.end() ) ++bb_atoms_found;
	if ( xyz.find( " C  " ) != xyz.end() ) ++bb_atoms_found;
	if ( xyz.find( " Pbb" ) != xyz.end() ) ++bb_atoms_found;

	if ( bb_atoms_found >= 1 && bb_atoms_found < 3 ) {
		is_l_aa = true;
		return;
	}
	// Positive angles are D
	core::Real characteristic_angle = 0;

	// Explicitly exclude peptoids and PNAs.
	if ( xyz.find( " CA " ) != xyz.end() && xyz.find( " CA1" ) == xyz.end() && xyz.find( " NG " ) == xyz.end() ) {
		// There are four atoms bonded to CA.
		if ( xyz.find( " Pbb" ) != xyz.end() ) {
			// Phosphonate
			characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " Pbb" ), xyz.at( " CB " ), xyz.at( " CA " ) );
		} else if ( xyz.find( " CM " ) != xyz.end() && name3 != "MLZ" ) { // methyllysine also uses CM, illustrating the weakness of this method.
			// beta
			if ( xyz.find( " CB " ) != xyz.end() ) {
				characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " CM " ), xyz.at( " CB " ), xyz.at( " CA " ) );
			} else if ( xyz.find( " CB1" ) != xyz.end() && xyz.find( " CB2" ) != xyz.end() ) {
				characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " CM " ), xyz.at( " CB1" ), xyz.at( " CB2" ) );
			} // other possibilities: B3G
		} else {
			// alpha
			if ( xyz.find( " CB " ) != xyz.end() ) {
				characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " C  " ), xyz.at( " CB " ), xyz.at( " CA " ) );
			} else if ( xyz.find( " CB1" ) != xyz.end() && xyz.find( " CB2" ) != xyz.end() ) {
				// CB1 is designated the L configuration controller.
				characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " C  " ), xyz.at( " CB1" ), xyz.at( " CA " ) );
			} // other possibilities: GLY
		}
		( characteristic_angle > 0 ) ? is_d_aa = true : is_l_aa = true;

		// Gammas--notably we need all this because just C2 C3 C4 are also had by sugars.
		// What a terrible method.
	} else if ( xyz.find( " C2 " ) != xyz.end() && xyz.find( " C3 " ) != xyz.end() && xyz.find( " C4 " ) != xyz.end() && xyz.find( " C  " ) != xyz.end() && xyz.find( " O  " ) != xyz.end() && xyz.find( " N  " ) != xyz.end() ) {
		// If we have a gamma, assign based on the stereo of the first carbon from C
		// This is an if else if NOT because we expect these to be mutually exclusive
		// but because CB3 only matters if there is no CB2.
		// We are NOT handling disubstituted.
		// AMW: Vikram's params maintain the position of the 2 in the name...
		if ( xyz.find( "CB2 " ) != xyz.end() ) {
			characteristic_angle = numeric::dihedral_degrees( xyz.at( " C3 " ), xyz.at( " C  " ), xyz.at( "CB2 " ), xyz.at( " C2 " ) );
		} else if ( xyz.find( "CB3 " ) != xyz.end() ) {
			characteristic_angle = numeric::dihedral_degrees( xyz.at( " C4 " ), xyz.at( " C2 " ), xyz.at( "CB3 " ), xyz.at( " C3 " ) );
		} else if ( xyz.find( "CB4 " ) != xyz.end() ) {
			characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " C3 " ), xyz.at( "CB4 " ), xyz.at( " C4 " ) );
		}
		( characteristic_angle > 0 ) ? is_d_aa = true : is_l_aa = true;
	}

	//if ( characteristic_angle == 0 ) {
	// TR.Warning << "No chiral information detected for " << name3 << std::endl;
	//}
}

bool
heavy_atom_names_match( ResidueType const & first, ResidueType const & second ) {
	if ( first.nheavyatoms() != second.nheavyatoms() ) {
		return false;
	}

	std::set< std::string > names;
	for ( core::Size ii(1); ii <= first.nheavyatoms(); ++ii ) {
		//if( first.atom_name(ii) == "OXT" ) { continue; }
		std::string name( first.atom_name(ii) );
		names.insert( utility::strip_whitespace( name ) );
	}

	for ( core::Size jj(1); jj <= second.nheavyatoms(); ++jj ) {
		std::string name( second.atom_name(jj) );
		if ( ! names.count( utility::strip_whitespace( name ) ) ) {
			//if( second.atom_name(jj) == "OXT" ) { continue; }
			TR.Debug << "Name mismatch: '" << second.atom_name(jj) <<"'" << std::endl;
			return false;
		}
	}

	// Since they have the same number of heavy atoms,
	// first can't have any names that second doesn't - at least without
	// second having one that first doesn't (assuming unique names)
	return true;
}


// Are these main-chain torsions also ring torsions?
bool
is_mainchain_torsion_also_ring_torsion( ResidueType const & res_type, uint torsion_index )
{
	using namespace std;
	using namespace utility;

	if ( res_type.is_cyclic() ) {
		vector1< uint > const mainchain_atom_indices( res_type.mainchain_atoms() );
		if ( torsion_index < mainchain_atom_indices.size() ) {
			pair< uint, uint > const mainchain_bond_atom_indices(
					mainchain_atom_indices[ torsion_index ], mainchain_atom_indices[ torsion_index + 1 ] );
			Size const n_ring_torsions( res_type.n_nus() );
			for ( uint i( 1 ); i <= n_ring_torsions; ++i ) {
				vector1< uint> const ring_torsion_definition( res_type.nu_atoms()[ i ] );
				pair< uint, uint > const nu_bond_atom_indices(
						ring_torsion_definition[ 2 ], ring_torsion_definition[ 3 ] );
				if ( mainchain_bond_atom_indices == nu_bond_atom_indices ) {
					return true;
				}
			}
		}
	}
	return false;
}

}  // namespace chemical
}  // namespace core
