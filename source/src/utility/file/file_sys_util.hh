// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/file_sys_util.hh
/// @brief  Platform independent operations on files (except I/O)
/// @author David Kim (dekim@u.washington.edu)
/// @author Ion Yannopoulos (ion@rosettacommons.org)


#ifndef INCLUDED_utility_file_file_sys_util_hh
#define INCLUDED_utility_file_file_sys_util_hh

// Utility Headers
#include <utility/file/FileName.hh>

// C++ headers
#include <fstream>

#include <utility/vector1.hh>

namespace utility {
namespace file {


/// @brief Does file exist?
bool
file_exists( std::string const & path );

/// @brief Is the file a directory?
bool
is_directory( std::string const & path );

/// @brief Delete file
int
file_delete( std::string const & path );


/// @brief Extension of a file name
std::string
file_extension( std::string const & filename );


/// @brief Prefix of a file name
std::string
file_basename( std::string const & filename );


/// @brief File size
long
file_size( std::string const & filename );

/// @brief current working directory
std::string
cwd();

/// @brief Create a blank file if it doesn't already exist
bool
create_blank_file( std::string const & blank_file );

/// @brief Find an unused random tempfile name with a given prefix (which may include a directory)
std::string
create_temp_filename( std::string const & dir, std::string const & prefix );

/// @brief Create a directory if it doesn't already exist
bool
create_directory(
	std::string const & dir_path
);


/// @brief Create a directory and its parent directories if they doesn't already exist
bool
create_directory_recursive(
	std::string const & dir_path
);


/// @brief Try to open file for read a few times just in case it is locked (from BOINC LIB)
bool
trytry_ifstream_open(
	std::ifstream & ifstream_,
	std::string const & name,
	std::ios_base::openmode open_mode
);


/// @brief Try to open file for write a few times just in case it is locked (from BOINC LIB)
bool
trytry_ofstream_open(
	std::ofstream & ofstream_,
	std::string const & name,
	std::ios_base::openmode open_mode
);

int list_dir (std::string dir, utility::vector1<std::string> & files);

FileName combine_names(utility::vector1<std::string> file_name_strings);

} // namespace file
} // namespace utility


#endif // INCLUDED_utility_file_file_sys_util_HH
