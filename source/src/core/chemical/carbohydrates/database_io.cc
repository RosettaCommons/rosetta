// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/carbohydrates/database_io.cc
/// @brief   Database input/output function definitions for carbohydrate-specific data.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/chemical/carbohydrates/database_io.hh>

// Project header
#include <core/types.hh>

// Utility headers
#include <utility/io/util.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// Basic header
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ header
#include <sstream>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.chemical.carbohydrates.database_io" );


namespace core {
namespace chemical {
namespace carbohydrates {

// Input //////////////////////////////////////////////////////////////////////

// Try various combinations to locate the specific glycan sequence file being requested by the user.
/// @details  The default directory to search is: database/chemical/carbohydrates/common_glycans/\n
/// The default file extension is: .iupac
std::string
find_glycan_sequence_file( std::string filename )
{
	using namespace utility::io;

	std::string dir = "chemical/carbohydrates/common_glycans/";
	std::string ext =  ".iupac";
	return basic::database::find_database_path( dir, filename, ext );
}

// Read a single-line glycan sequence file.
std::string
read_glycan_sequence_file( std::string filename )
{
	utility::vector1< std::string > const lines( utility::io::get_lines_from_file_data( filename ) );
	if ( lines.size() != 1 ) {
		utility_exit_with_message( "A glycan sequence file must contain a single line of text." );
	}
	return lines.front();
}


// Check if a string variable's value is "N/A", and if so, sets it to "".
void
check_if_applicable( std::string & variable )
{
	if ( variable == "N/A" ) { variable = ""; }
}

void
replace_underscores_with_spaces( std::string & phrase )
{
	replace( phrase.begin(), phrase.end(), '_', ' ' );
}


// Return a map of strings to strings, which are saccharide-specific 3-letter codes mapped to IUPAC roots, read from a
// database file.
std::map< std::string, std::pair< std::string, char > >
read_codes_and_roots_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	map< string, pair< string, char > > codes_to_roots;

	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {
		istringstream line_word_by_word( lines[ i ] );
		string key;  // The map key is the 3-letter code, e.g., "Glc", for "glucose".
		string root;  // The IUPAC root, e.g., "gluc", for " glucose".
		char stereochem;  // The default stereochemistry, L or D, or * for inherent.

		line_word_by_word >> key >> root >> stereochem;

		codes_to_roots[ key ] = make_pair( root, stereochem );
	}

	TR.Debug << "Read " << codes_to_roots.size() << " 3-letter code mappings from the carbohydrate database." << endl;

	return codes_to_roots;
}

// Return a map of Sizes to pairs of char and string, which are ring sizes mapped to 1-letter affixes and morphemes,
// respectively, read from a database file.
std::map< core::Size, std::pair< char, std::string > >
read_ring_sizes_and_morphemes_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	map< Size, pair< char, string > > ring_size_to_morphemes;

	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {
		istringstream line_word_by_word( lines[ i ] );
		Size key;  // The map key is the ring size.
		char affix;  // The first element of the pair is a 1-letter affix, e.g., "f", for a "furanose".
		string morpheme;  // The second element of the pair is the internal morpheme, e.g., "ofuran", for a "furanose".

		line_word_by_word >> key >> affix >> morpheme;

		if ( key < 3 ) {
			utility_exit_with_message( "read_ring_sizes_and_morphemes_from_database_file: "
				"invalid ring size; rings cannot have less than 3 atoms!" );
		}
		if ( affix == 'X' ) { affix = '\0'; }  // Some ring sizes don't have accepted affixes.

		ring_size_to_morphemes[ key ] = make_pair( affix, morpheme );
	}

	TR.Debug << "Read " << ring_size_to_morphemes.size() <<
		" ring size mappings from the carbohydrate database." << endl;

	return ring_size_to_morphemes;
}

// Return a table of nomenclature data for sugar modifications, read from a database file.
SugarModificationsNomenclatureTable
read_nomenclature_table_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	SugarModificationsNomenclatureTable table;

	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {
		istringstream line_word_by_word( lines[ i ] );
		std::string key;  // The map key is the name of the type of modification, derived from the VariantType name.
		SugarModificationsNomenclatureTableRow row;

		line_word_by_word >> key >> row.substituent_full_name >> row.implies >> row.short_affix >> row.patch_name >>
			row.default_position >> row.has_inherent_position /*>> row.reducing_end_suffix >> glycoside_suffix >>
			row.following_word_or_phrase*/;

		replace_underscores_with_spaces( key );
		check_if_applicable( row.substituent_full_name );
		check_if_applicable( row.implies );
		check_if_applicable( row.short_affix );
		check_if_applicable( row.patch_name );
		//check_if_applicable( row.reducing_end_suffix );
		//check_if_applicable( row.glycoside_suffix );
		//check_if_applicable( row.following_word_or_phrase );

		// TODO: After ResidueProperties static const data is moved into a singleton, check that each key is in fact
		// derived from accepted ResiduePropertys

		table[ key ] = row;
	}

	TR.Debug << "Read " << table.size() <<
		" rows from the sugar modifications table in the carbohydrate database." << endl;

	return table;
}

// Return a map of linkage conformer data, read from a database file.
/// @author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)
LinkageConformers
read_linkage_conformers_from_database_file( std::string const & filename ) {
	using namespace std;
	using namespace utility;
	using namespace core::id;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	LinkageConformers conformer_data_structure;


	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {
		istringstream line_word_by_word( lines[ i ] );
		TR.Debug << lines[ i ] << std::endl;

		// Every line is a new conformer.
		LinkageConformerData conformer;
		string non_red_res, red_res;
		string conformer_bins;
		Real pop;
		Angle phi_mean, phi_sd, psi_mean, psi_sd, omega_mean, omega_sd;

		line_word_by_word >> non_red_res >> red_res >> conformer_bins >> pop >> phi_mean >> phi_sd >> psi_mean >> psi_sd;

		conformer.population = pop;
		conformer.mean_sd.push_back( make_pair( phi_mean, phi_sd ) );
		conformer.mean_sd.push_back( make_pair( psi_mean, psi_sd ) );
		conformer.conformer_bins = conformer_bins;

		utility::vector1< std::string > entries = utility::string_split_multi_delim(lines[i]);
		core::Size total_entries = entries.size();
		core::Size total_omegas = (total_entries - 8) / 2;

		TR.Debug << "split " <<utility::to_string( entries ) << std::endl;
		TR.Debug << "total entries " << total_entries << std::endl;
		TR.Debug << "total omegas" << total_omegas << std::endl;

		if ( total_omegas > 0 ) {
			for ( core::Size omega_n = 1; omega_n <= total_omegas; ++ omega_n ) {

				line_word_by_word >> omega_mean;
				line_word_by_word >> omega_sd;
				if ( line_word_by_word.fail() ) {
					utility_exit_with_message( "read_linkage_conformers_from_database_file: "
						"Omega angles mush be listed in mean/sd pairs!" );
				}
				conformer.omega_mean_sd.push_back( make_pair( omega_mean, omega_sd ) );
				TR.Debug << "Omega "<< omega_n << omega_sd << std::endl;
			}
		}

		pair< string, string > const linkage_pair( make_pair( non_red_res, red_res ) );
		conformer_data_structure[ linkage_pair ].push_back( conformer );
	}

	return conformer_data_structure;
}

std::map< std::string, std::string >
read_short_names_to_iupac_format_string( std::string const & dir, std::string common_mapping_path)
{
	using utility::vector1;
	using std::string;

	std::map< string, string > short_names_to_iupac_string;

	vector1< string > const lines( utility::io::get_lines_from_file_data( common_mapping_path ) );

	for ( Size i = 1; i <= lines.size(); ++i ) {
		string line = lines[ i ];

		TR << line << std::endl; //Debugging for now.


		std::istringstream line_word_by_word( line );

		std::string short_name;
		std::string iupac_file_name;
		line_word_by_word >> short_name >> iupac_file_name;

		std::string iupac_path = basic::database::find_database_path( dir, iupac_file_name);
		std::string iupac_sequence = read_glycan_sequence_file( iupac_path );

		short_names_to_iupac_string[ short_name ] = iupac_sequence;

	}

	return short_names_to_iupac_string;


}


}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
