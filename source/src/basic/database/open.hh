// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/database/open.hh
/// @brief  Functions for opening database files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_basic_database_open_hh
#define INCLUDED_basic_database_open_hh


// Utility headers
#include <utility/io/izstream.fwd.hh>

// C++ headers
#ifdef WIN32
#include <string>
#else
#include <iosfwd>
#endif

namespace basic {
namespace database {

/// @brief Does a database file exist?
/* Undefinded, commenting out to fix PyRosetta build
bool
exists(
std::string const & db_file
); */

/// @brief Open a database file on a provided stream
bool
open(
	utility::io::izstream & db_stream,
	std::string const & db_file,
	bool warn = true
);

/// @brief Full-path database file name
std::string
full_name(
	std::string const & db_file,
	bool warn = true
);

/// @brief Does cache file (absolute path) exist?
/// if dir_only is true, will return true if the cache file could be created.
bool
find_cache_file(
	std::string const & cache_file,
	bool dir_only
);

/// @brief Get the (absolute) path to a given cached file.
/// If source_file is given, it's the full path to the source database file that's being cached.
/// If for_writing is true, will only check that the given file would be creatable.
/// Will return an empty string if it can't find a cache file.
std::string
full_cache_name(
	std::string const & short_name,
	std::string const & source_file,
	bool for_writing
);

} // namespace database
} // namespace basic


#endif // INCLUDED_basic_io_database_open_HH
