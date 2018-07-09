// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/NomenclatureManager.hh
/// @brief   Method definitions for NomenclatureManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/NomenclatureManager.hh>

// Package header
#include <core/io/alt_codes_io.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/io/izstream.hh>
#include <utility/io/util.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>
#include <boost/algorithm/string.hpp>


// Construct tracer.
static basic::Tracer TR( "core.io.pdb.NomenclatureManager" );


namespace core {
namespace io {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
/// @return  If a file of alternative 3-letter codes has been provided at the command line, return a pair of Rosetta
/// names (3-letter code and base ResidueType name, if available) from the given PDB 3-letter code.  If no file has
/// been provided, or if there is no alternative for the given code, simply return the PDB code and an empty string.
std::pair< std::string, std::string >
NomenclatureManager::rosetta_names_from_pdb_code( std::string const & pdb_code )
{
	using namespace std;
	using namespace basic::options;

	// We used to check for option use here -- but there is now a default use.
	// This is because we need -- for ligand_water_docking -- to recognize
	// both TP3 and WAT as waters.
	AltCodeMap const & alt_codes( get_instance()->get_alternate_3_letter_code_map() );
	AltCodeMap::const_iterator alt_code_key_value_pair( alt_codes.find( pdb_code ) );
	if ( alt_code_key_value_pair == alt_codes.end() ) {  // Is there an alternate for this code?
		return make_pair( pdb_code, "" );  // If there is no alternate, simply return what was passed.
	}

	// Get the value (a tuple) of this key/value pair and tie it to temp. variables.
	string const & name3( get< 0 >( alt_code_key_value_pair->second ) );
	string const & base_name( get< 1 >( alt_code_key_value_pair->second ) );
	// We use another function, default_mainchain_connectivity_from_pdb_code, to access the 3rd (index 2) value.
	// We do not need the 4th (index 3) value for input purposes, only output.

	TR.Debug << "Accepting alternate code " << pdb_code << " for " << name3 << '.' << endl;
	return make_pair( name3, base_name );
}

/// @return  If a file of alternative 3-letter codes has been provided at the command line, return the default main-
/// chain connectivity (if available) from the given PDB 3-letter code.  If no file has been provided, or if there is
/// no default value for the given code, simply return a null character.
char
NomenclatureManager::default_mainchain_connectivity_from_pdb_code( std::string const & pdb_code )
{
	using namespace basic::options;

	AltCodeMap const & alt_codes( get_instance()->get_alternate_3_letter_code_map() );
	AltCodeMap::const_iterator alt_code_key_value_pair( alt_codes.find( pdb_code ) );
	if ( alt_code_key_value_pair == alt_codes.end() ) {  // Is there an alternate for this code?
		return '\0';    // If there is no alternate, simply return a null char.
	}

	// Get the value (a tuple) of this key/value pair and grab the connectivity (3rd value, index 2).
	return std::get< 2 >( alt_code_key_value_pair->second );
}

// This function returns the pdb code that corresponds to the residue name.
/// @return  The unique PDB 3-letter code that corresponds to the give ResidueType name or an empty string if no match
/// is found.
/// @author  Brandon Frenz
/// @author  Labonte <JWLabonte@jhu.edu>
std::string
NomenclatureManager::pdb_code_from_rosetta_name( std::string const & rosetta_name )
{
	using namespace std;
	using namespace utility;

	debug_assert( rosetta_name != "" );

	// TODO: We need a better way of getting these patch names, as hypothetically, other PTMs will have the same issue.
	// The list below is exclusively for saccharide ResidueTypes, and this function should work with non-saccharides.
	vector1< string > const non_modification_patches = {
		"non-reducing_end",
		"reducing_end",
		"cutpoint_lower",
		"cutpoint_upper",
		"->1)-branch",
		"->2)-branch",
		"->3)-branch",
		"->4)-branch",
		"->5)-branch",
		"->6)-branch",
		"->7)-branch",
		"->8)-branch",
		"->9)-branch"
		};

	TR.Debug << "Looking for PDB code for " << rosetta_name << std::endl;

	// Separate the base name from the patches and main-chain designation, if necessary.
	vector1< string > base_name_and_patches;
	vector1< string > patches;
	if ( rosetta_name.find( ":" ) != string::npos ) {
		boost::split( base_name_and_patches, rosetta_name, []( char c ) { return c == ':'; } );
	} else {
		base_name_and_patches.push_back( rosetta_name );
	}

	// Start at 2 since #1 is the base name.
	for ( core::Size i = 2; i <= base_name_and_patches.size(); ++i ) {
		if ( ! non_modification_patches.has_value( base_name_and_patches[ i ] ) ) {
			patches.push_back( base_name_and_patches[ i ] );
		}
	}
	string base_name( base_name_and_patches[ 1 ] );
	if ( base_name.find( "->" ) != std::string::npos ) {
		// Remove the main-chain designation, if present, from the base name.
		vector1< string > base_name_vect;
		boost::split( base_name_vect, base_name, []( char c ) { return c == ')'; } );
		base_name = base_name_vect[ 2 ];  // e.g., -alpha-D-Glcp
	}

	// Loop over the alt code map.
	AltCodeMap const & alt_codes( get_instance()->get_alternate_3_letter_code_map() );
	for ( auto const & alt_code : alt_codes ) {
		// Find the base name of this AltCodeMap record.
		string alt_code_base_name( get< 1 >( alt_code.second ) );
		if ( alt_code_base_name.find( "->" ) != std::string::npos ) {
			vector1< string > alt_code_base_name_vect;
			boost::split( alt_code_base_name_vect, alt_code_base_name, []( char c ) { return c == ')'; } );
			alt_code_base_name = alt_code_base_name_vect[ 2 ];  // e.g., -alpha-D-Glcp
		}

		TR.Trace << " Base name for PDB code " << alt_code.first << ": " << alt_code_base_name << endl;

		if ( base_name == alt_code_base_name ) {
			TR.Debug << " Base names match!" << std::endl;

			vector1< string > const & alt_code_patches( get< 3 >( alt_code.second ) );

			if ( TR.Debug.visible() ) {
				TR.Debug << "  Patches for PDB code " << alt_code.first << ": ";
				for ( auto const & patch : alt_code_patches ) {
					TR.Debug << patch << ' ';
				}
				TR.Debug << std::endl;
			}

			if ( patches.size() == alt_code_patches.size() ) {
				bool has_all( true );
				for ( auto const & patch : alt_code_patches ) {
					if ( ! patches.has_value( patch ) ) {
						has_all = false;
						break;
					}
				}
				if ( has_all ) {
					TR.Debug << " Perfect match!" << endl;
					return alt_code.first;
				}
			}
		}
	}
	return "";  // If we get here, no match was found.
}

std::string
NomenclatureManager::annotated_sequence_from_modomics_oneletter_sequence( std::string const & seq )  {
	std::string annotated_seq;
	for ( char const c : seq ) {
		// Don't translate ACGU in the annotated_seq style
		if ( c == 'A' ) {
			annotated_seq += "a";
		} else if ( c == 'C' ) {
			annotated_seq += "c";
		} else if ( c == 'G' ) {
			annotated_seq += "g";
		} else if ( c == 'U' ) {
			annotated_seq += "u";
		} else if ( c == '#' ) {
			// For some reason, this isn't making it into the map -- despite code to ensure that it will
			// by only taking digraphs as commnets.
			annotated_seq += "X[OMG]";
		} else {
			if ( get_instance()->modomics_map().find( c ) == get_instance()->modomics_map().end() ) {
				TR.Error << "In sequence \"" << seq << "\", character \'" << c << "\' was not recognized." << std::endl;
				utility_exit_with_message( "A character in your provided modomics-style sequence was not recognized." );
			} else {
				annotated_seq += "X[" + get_instance()->modomics_map().at( c ) + "]";
			}
		}
	}
	return annotated_seq;
}

std::string
NomenclatureManager::annotated_sequence_from_IUPAC_sequence( std::string const & seq ) {
	std::string annotated_seq;
	utility::vector1< std::string > parts = utility::string_split( seq, '_' );
	for ( auto const & part : parts ) {
		if ( get_instance()->iupac_map().find( part ) == get_instance()->iupac_map().end() ) {
			TR.Error << "In sequence \"" << seq << "\", component\'" << part << "\' was not recognized." << std::endl;
			utility_exit_with_message( "A character in your provided IUPAC-style sequence was not recognized." );
		}
		annotated_seq += "X[" + get_instance()->iupac_map().at( part ) + "]";
	}
	return annotated_seq;
}


bool NomenclatureManager::is_NA( std::string const & name3 ) {
	std::set< std::string > const & na_set( get_instance()->na_set() );
	return na_set.count( name3 );
}

bool NomenclatureManager::is_old_RNA( std::string const & name3 ) {
	std::set< std::string > const & old_rna_set( get_instance()->old_rna_set() );
	return old_rna_set.count( name3 );
}

bool NomenclatureManager::is_old_DNA( std::string const & name3 ) {
	std::set< std::string > const & old_dna_set( get_instance()->old_dna_set() );
	return old_dna_set.count( name3 );
}

bool NomenclatureManager::decide_is_d_aa( std::string const & name3 ) {
	std::set< std::string > const & d_aa_set( get_instance()->d_aa_set() );
	return d_aa_set.count( name3 );
}

bool NomenclatureManager::decide_is_l_aa( std::string const & name3 ) {
	std::set< std::string > const & l_aa_set( get_instance()->l_aa_set() );
	return l_aa_set.count( name3 );
}

bool NomenclatureManager::decide_is_known_achiral( std::string const & name3 ) {
	std::set< std::string > const & achiral_set( get_instance()->achiral_set() );
	return achiral_set.count( name3 );
}

bool NomenclatureManager::is_metal( std::string const & name3 ) {
	std::set< std::string > const & metal_set( get_instance()->metal_set() );
	return metal_set.count( name3 );
}

bool NomenclatureManager::is_sugar( std::string const & name3 ) {
	std::set< std::string > const & sugar_set( get_instance()->sugar_set() );
	return sugar_set.count( name3 );
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
NomenclatureManager::NomenclatureManager()
{
	using namespace std;
	using namespace basic::options;
	using namespace core;

	utility::vector1< string > const & file_list( option[ OptionKeys::in::alternate_3_letter_codes ]() );
	Size const n_files( file_list.size() );
	for ( uint i( 1 ); i <= n_files; ++i ) {
		string const & filename( find_alternate_codes_file( file_list[ i ] ) );
		AltCodeMap alt_codes( read_alternative_3_letter_codes_from_database_file( filename ) );
		alt_codes_.insert( alt_codes.begin(), alt_codes.end() );
	}

	utility::vector1< string > const & file_list2 = option[ OptionKeys::in::name3_property_codes ]();
	for ( Size jj = 1; jj <= file_list2.size(); ++jj ) {
		utility::vector1< string > const lines( utility::io::get_lines_from_file_data( find_alternate_codes_file( file_list2[ jj ] ) ) );
		for ( Size ii = 1; ii <= lines.size(); ++ii ) {

			if ( lines[ ii ].size() == 0 || lines[ ii ][ 0 ] == '#' ) continue;
			istringstream word_by_word( lines[ ii ] );

			string name3, value;
			getline( word_by_word, value, '\t' );
			getline( word_by_word, name3, '\t' );

			if ( value == "NA" ) {
				is_NA_.insert( name3 );
			} else if ( value == "OLD_DNA" ) {
				is_old_DNA_.insert( name3 );
			} else if ( value == "OLD_RNA" ) {
				is_old_RNA_.insert( name3 );
			} else if ( value == "L_AA" ) {
				l_aa_set_.insert( name3 );
			} else if ( value == "D_AA" ) {
				d_aa_set_.insert( name3 );
			} else if ( value == "ACHIRAL" ) {
				achiral_set_.insert( name3 );
			} else if ( value == "METAL" ) {
				metal_set_.insert( name3 );
			} else if ( value == "SUGAR" ) {
				sugar_set_.insert( name3 );
			} else {
				utility_exit_with_message( "Line in name3 properties file \"" + file_list2[ jj ] + "\" malformed:\"" + lines[ ii ] + "\"" );
			}
		}
	}

	// Set up annotated sequence generation from:
	// Modomics 1-letter codes
	// IUPAC designations
	//   NOTE: We have to use the #/ digraph to indicate comments because # is
	//   one of the Modomics symbols.
	//   NOTE: Moving to a single file for both of these purposes because we don't
	//   want to have to update two places with annotated sequence entries
	utility::vector1< string > const iupac_lines( utility::io::get_lines_from_file_data( basic::database::full_name( "input_output/modomics_and_iupac_to_ann.txt" ) ) );
	for ( auto const line : iupac_lines ) {
		if ( line.size() == 0 || ( line[ 0 ] == '#' && line[ 1 ] == '/' ) ) continue;
		istringstream word_by_word( line );
		std::string modomics, iupac, ann;
		std::getline( word_by_word, modomics, '\t' );
		std::getline( word_by_word, iupac, '\t' );
		std::getline( word_by_word, ann, '\t' );

		// Only first char of modomics string (well, it's a single char string)
		annotated_seq_from_modomics_map_[ modomics[0] ] = ann;
		annotated_seq_from_IUPAC_map_[ iupac ] = ann;
	}
}

// Get the map requested, creating it if necessary.
// Called by the public static method rosetta_names_from_pdb_code()
AltCodeMap const &
NomenclatureManager::get_alternate_3_letter_code_map() const
{
	return alt_codes_;
}

// Try various combinations to locate the specific file being requested by the user.
// (inspired by core::scoring::ScoreFunction::find_weights_file())
std::string
NomenclatureManager::find_alternate_codes_file( std::string const & filename )
{
	using namespace utility::io;

	std::string const & path( basic::database::full_name( "input_output/3-letter_codes/" ) );
	std::string const ext( ".codes" );

	izstream potential_file( filename );
	if ( potential_file.good() ) {
		return filename;
	} else {
		izstream potential_file( filename + ext );  // Perhaps the user didn't use the .codes extension.
		if ( potential_file.good() ) {
			return filename + ext;
		} else {
			izstream potential_file( path + filename);  // Let's assume it's in the database in the usual spot.
			if ( potential_file.good() ) {
				return path + filename;
			} else {
				izstream potential_file( path + filename + ext );  // last try
				if ( potential_file.good() ) {
					return path + filename + ext;
				} else {
					utility_exit_with_message( "Unable to open alternative codes file. Neither ./" + filename +
						" nor " + "./" + filename + ext +
						" nor " + path + filename +
						" nor " + path + filename + ext + " exists." );
				}
			}
		}
	}
	return "INCONCEIVABLE!";  // Code can never reach here.
}


}  // namespace io
}  // namespace core
