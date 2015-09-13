// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

// Basic header
#include <basic/Tracer.hh>

// C++ header
#include <sstream>


// Construct tracer.
static thread_local basic::Tracer TR( "core.chemical.carbohydrates.database_io" );


namespace core {
namespace chemical {
namespace carbohydrates {

// Return a map of strings to strings, which are saccharide-specific 3-letter codes mapped to IUPAC roots, read from a
// database file.
std::map< std::string, std::string >
read_codes_and_roots_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	map< string, string > codes_to_roots;

	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {
		istringstream line_word_by_word( lines[ i ] );
		string key;  // The map key is the 3-letter code, e.g., "Glc", for "glucose".
		string value;  // The map value is the IUPAC root, e.g., "gluc", for " glucose".

		line_word_by_word >> key >> value;

		codes_to_roots[ key ] = value;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << codes_to_roots.size() << " 3-letter code mappings from the carbohydrate database." << endl;
	}

	return codes_to_roots;
}

// Return a map of Sizes to pairs of char and string, which are ring sizes mapped to 1-letter affixes and morphemes,
// respectively.
std::map< core::Size, std::pair< char, std::string > >
read_ring_sizes_and_morphemes_fromt_database_file( std::string const & filename )
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
			utility_exit_with_message( "read_ring_sizes_and_morphemes_fromt_database_file: "
				"invalid ring size; rings cannot have less than 3 atoms!" );
		}
		if ( affix == 'X' ) { affix = '\0'; }  // Some ring sizes don't have accepted affixes.

		ring_size_to_morphemes[ key ] = make_pair( affix, morpheme );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << ring_size_to_morphemes.size() <<
			" ring size mappings from the carbohydrate database." << endl;
	}

	return ring_size_to_morphemes;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
