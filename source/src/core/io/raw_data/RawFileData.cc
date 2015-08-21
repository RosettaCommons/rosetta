// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/RawFileData.cc
///
/// @brief FileData base class
/// @author James Thompson, Monica Berrondo

// C++ Headers
#include <vector>
#include <string>
#include <map>

// mini headers
#include <core/io/raw_data/Raw.fwd.hh>
#include <core/io/raw_data/RawFileData.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>


namespace core {
namespace io {
namespace raw_data {

utility::vector1< std::string > RawFileData::read_tags_fast( std::string const filename ) const {
	utility::vector1< std::string > tags_in_file;
	utility::io::izstream data( filename.c_str() );
	if ( !data ) {
		std::cerr << "ERROR:: Unable to open silent_input file: " << filename << std::endl;
		return tags_in_file;
	}

	std::string line;
	getline( data, line ); // sequence line
	getline( data, line ); // score line

	while ( getline(data,line) ) {
		if ( line.substr(0,7) == "SCORE: " ) {
			std::istringstream l( line );

			std::string current_word;
			while ( !l.fail() ) {
				l >> current_word;
			}
			tags_in_file.push_back( current_word );
		}
	} // while( getline(data,line) )

	return tags_in_file;
} // read_tags_fast

void RawFileData::write_all(
	const std::string filename,
	std::map < std::string, core::Real > const & score_map
) {
	bool print_header( true );
	utility::io::ozstream output;
	if ( !utility::file::file_exists( filename ) ) {
		output.open( filename );
		print_header = true;
	} else {
		output.open_append( filename );
	}

	for ( core::io::raw_data::RawFileData::iterator iter = begin(), it_end = end(); iter != it_end; ++iter ) {
		if ( print_header ) {
			iter->print_header( output, score_map );
			print_header = false;
		}

		iter->print_scores      ( output, score_map );
		iter->print_conformation( output );
	}
} // write_all

} // namespace silent
} // namespace io
} // namespace core
