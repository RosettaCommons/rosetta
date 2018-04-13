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
/// names (3-letter code and base ResidueType name, if available); if no file has been provided, or if there is no
/// alternative for the given code, simply return the PDB code and an empty string.
std::pair< std::string, std::string >
NomenclatureManager::rosetta_names_from_pdb_code( std::string const & pdb_code )
{
	using namespace std;
	using namespace basic::options;

	if ( option[ OptionKeys::in::alternate_3_letter_codes ].active() ) {  // Are alternate codes allowed?
		AltCodeMap const & alt_codes( get_instance()->get_alternate_3_letter_code_map() );
		auto alt_code_tuple( alt_codes.find( pdb_code ) );
		if ( alt_code_tuple != alt_codes.end() ) {  // Is there an alternate for this code?
			// Get the value of this key/value pair. This function now converts to a pair from the tuple in order to maintain backward compatability. BF.
			std::string name3 = std::get<0>(alt_code_tuple->second);
			std::string base_name = std::get<1>(alt_code_tuple->second);
			pair< string, string > const & rosetta_names = std::make_pair( name3, base_name );
			TR.Debug << "Accepting alternate code " << pdb_code << " for " << rosetta_names.first << '.' << endl;
			return rosetta_names;
		}
	}
	return std::make_pair( pdb_code, "" );
}

//This function returns the pdb code that corresponds to the residue name. It will fail if no match can be found. If this happens add the patches
//you want to match to the pdb_codes database.
std::string
NomenclatureManager::pdb_code_from_rosetta_name(std::string const & rosetta_name ){
	utility::vector1<std::string> split_vect;
	boost::split(split_vect, rosetta_name, [](char c){return c == ':';});
	utility::vector1<std::string> patches;
	//start at 2 since #1 is the base name
	for ( core::Size i=2; i<=split_vect.size(); i++ ) {
		patches.push_back(split_vect[i]);
	}
	std::string base_name;
	utility::vector1<std::string> base_name_vect;
	boost::split(base_name_vect, split_vect[1], [](char c){ return c == ')';});
	base_name = base_name_vect[2];
	//Add all patches that don't code for sugars here
	utility::vector1<std::string> base_patches;
	base_patches.push_back("non-reducing_end");
	base_patches.push_back("reducing_end");
	base_patches.push_back("branch_lower_terminus");
	base_patches.push_back("cutpoint_lower");
	base_patches.push_back("cutpoint_upper");
	for ( core::Size i=1; i<=patches.size(); i++ ) {
		if ( patches[i].find("branch") != std::string::npos ) {
			base_patches.push_back(patches[i]);
		}
	}

	//loop over the alt code map
	AltCodeMap const & alt_codes( get_instance()->get_alternate_3_letter_code_map() );
	core::Size max_patches_matched = 0;
	std::string best_match = "";
	utility::vector1<std::string> best_patch_match;
	for ( auto const & alt_code : alt_codes ) {
		//determine whether the linkage number is related to the name
		bool match_all = false;
		if ( std::get<1>(alt_code.second).find("->?") != std::string::npos ) {
			match_all = true;
		}
		utility::vector1<std::string> matching_patches = std::get<2>(alt_code.second);

		//find the base name
		std::string base_tag;
		if ( match_all ) {
			utility::vector1<std::string> base_tag_vect;
			boost::split(base_tag_vect, std::get<1>(alt_code.second), [](char c){ return c == ')';});
			base_tag = base_tag_vect[2];
		} else {
			base_tag = split_vect[1];
		}

		if ( base_name == base_tag ) {
			bool has_all = true;
			core::Size match_count = 0;
			for ( core::Size i=1; i<=matching_patches.size(); i++ ) {
				if ( matching_patches[i] == "default" ) continue;
				if ( !patches.has_value(matching_patches[i]) ) {
					has_all = false;
				} else {
					match_count+=1;
				}
			}
			if ( has_all == false ) break;
			if ( match_count > max_patches_matched ) {
				max_patches_matched = match_count;
				best_match = alt_code.first;
				best_patch_match = matching_patches;
			}
			if ( best_match == "" && matching_patches.has_value("default") ) {
				best_match = alt_code.first;
			}
		}
	}
	//kill the dump if any patches aren't recognized
	bool no_match = false;
	for ( core::Size i=1; i<=patches.size(); i++ ) {
		if ( !best_patch_match.has_value(patches[i]) && !base_patches.has_value(patches[i]) ) {
			no_match = true;
		}
	}
	if ( best_match == "" ) {
		no_match = true;
	}
	if ( no_match ) {
		TR.Warning << "Could not match " << rosetta_name << " to pdb code. Please turn off the write_glycan_pdb_code flag. Developers should add the full patch name to the pdb codes database" << std::endl;
		runtime_assert(no_match == false);
	}
	return best_match;
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

	if ( option[ OptionKeys::in::alternate_3_letter_codes ].user() ) {
		utility::vector1< string > const & file_list( option[ OptionKeys::in::alternate_3_letter_codes ]() );
		Size const n_files( file_list.size() );
		for ( uint i( 1 ); i <= n_files; ++i ) {
			string const & filename( find_alternate_codes_file( file_list[ i ] ) );
			AltCodeMap alt_codes( read_alternative_3_letter_codes_from_database_file( filename ) );
			alt_codes_.insert( alt_codes.begin(), alt_codes.end() );
		}
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
