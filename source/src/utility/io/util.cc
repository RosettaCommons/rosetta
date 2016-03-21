// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    utility/io/util.cc
/// @brief   General database input/output utility function definitions.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <utility/io/util.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>

namespace utility {
namespace io {

// General method that opens a file and returns its data as a list of lines after checking for errors.
/// @details  Blank and commented lines are not returned and the file is closed before returning the lines.
/// @author   Labonte <JWLabonte@jhu.edu>
utility::vector1< std::string >
get_lines_from_file_data( std::string const & filename )
{
	using namespace std;
	using namespace utility;
	using namespace utility::file;
	using namespace utility::io;

	// Check if file exists.
	if ( ! file_exists( filename ) ) {
		utility_exit_with_message( "Cannot find database file: '" + filename + "'" );
	}

	// Open file.
	izstream data( ( filename.c_str() ) );
	if ( ! data.good() ) {
		utility_exit_with_message( "Unable to open database file: '" + filename + "'" );
	}

	string line;
	vector1< string > lines;

	while ( getline( data, line ) ) {
		
		trim( line, " \t\n" );  // Remove leading and trailing whitespace.
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) { continue; }  // Skip comments and blank lines.
		lines.push_back( line );
	}

	data.close();

	return lines;
}

}  // namespace io
}  // namespace utility
