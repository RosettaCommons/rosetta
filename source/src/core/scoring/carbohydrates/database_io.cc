// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/carbohydrates/database_io.cc
/// @brief   Database input/output function definitions for carbohydrate-specific scoring data.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/scoring/carbohydrates/database_io.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/io/util.hh>

// Basic header
#include <basic/Tracer.hh>

// C++ header
#include <sstream>


// Construct tracer.
static thread_local basic::Tracer TR( "core.scoring.carbohydrates.database_io" );


namespace core {
namespace scoring {
namespace carbohydrates {

// Return a table of Gaussian parameters read from a database file.
std::map< char, utility::vector1< Real > >
read_Gaussian_parameters_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	map< char, vector1< Real > > parameters;

	Size const n_lines( lines.size() );
	for ( uint line_num( 1 ); line_num <= n_lines; ++line_num ) {
		istringstream line_word_by_word( lines[ line_num ] );
		char name;  // The map key is the one-letter variable name of the parameter.
		Real value;
		vector1< Real > values;

		line_word_by_word >> name >> value;
		while ( ! line_word_by_word.fail() ) {
			values.push_back( value );
			line_word_by_word >> value;
		}

		parameters[ name ] = values;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << parameters.size() <<
			" Gaussian parameter sets from the carbohydrate scoring database." << endl;
	}

	return parameters;
}

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core
