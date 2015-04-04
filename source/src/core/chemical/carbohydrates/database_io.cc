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

// Utility header
#include <utility/io/util.hh>

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

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
