// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/gzip_util.hh
/// @brief  gzip utility functions
/// @author David Kim (dekim@u.washington.edu)


#ifndef INCLUDED_utility_file_gzip_util_hh
#define INCLUDED_utility_file_gzip_util_hh


// C++ headers
// Issue with Windows PyRosetta requires full header, no ifdef
#include <string>


namespace utility {
namespace file {


/// @brief gzip: file compression
long
gzip(
	std::string const & uncompressedfile,
	bool overwrite = false
);


/// @brief gunzip: file decompression
long
gunzip(
	std::string const & compressedfile,
	bool overwrite = false
);


} // namespace file
} // namespace utility


#endif // INCLUDED_utility_file_gzip_util_HH
