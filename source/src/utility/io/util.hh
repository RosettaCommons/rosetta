// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/util.hh
/// @brief  General database input/output utility functions


#ifndef INCLUDED_utility_io_util_hh
#define INCLUDED_utility_io_util_hh

#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

namespace utility {
namespace io {

template< class T >
void
read_vector( std::istream & is, vector1< T > & vec )
{
	vec.clear();
	T val;
	while ( is >> val ) {
		vec.push_back( val );
	}
}

template< class T >
void
write_vector( std::ostream & out, vector1< T > const & vec )
{
	for ( typename vector1< T >::const_iterator it = vec.begin(), eit = vec.end(); it != eit; ++it ) {
		out << *it << "\n";
	}
}

template< class T >
void
write_vector( std::string filename, vector1< T > const & vec )
{
	utility::io::ozstream out( filename );
	write_vector( out, vec );
}

/// @brief  General method that opens a file and returns its data as a list of lines after checking for errors.
utility::vector1< std::string >
get_lines_from_file_data( std::string const & filename );


}
}

#endif
