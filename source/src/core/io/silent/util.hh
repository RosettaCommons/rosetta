// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/util.hh
///
/// @brief utility functions for silent-files
/// @author James Thompson
/// @author Rhiju Das

#ifndef INCLUDED_core_io_silent_util_hh
#define INCLUDED_core_io_silent_util_hh

#include <map>
#include <string>
#include <utility/vector1.fwd.hh>

namespace core {
namespace io {
namespace silent {

/// @brief gzip all of the files in -out::file::silent().
void
gzip( void );

/// @brief figures out which decoys are already stored in silent file.
std::map< std::string, bool >
initialize_tag_is_done( std::string const & silent_file );

/// @brief needed for reading out RESNUM lines from silent files
void
figure_out_residue_numbers_from_line( std::istream & line_stream,
	utility::vector1< int > & residue_numbers,
	utility::vector1< char > & chains,
	utility::vector1< std::string > & segids );

/// @brief Changes blah.out to blah_LORES.out (if tag is "_LORES")
std::string
get_outfile_name_with_tag( std::string const & silent_file, std::string const & tag );

/// @brief used with -overwrite flag
void
remove_silent_file_if_it_exists( std::string const & silent_file);

} // namespace silent
} // namespace io
} // namespace core

#endif
