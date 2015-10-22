// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/boinc_util.hh
/// @author David Kim (dekim@u.washington.edu)

// Unit header
#include <utility/boinc/boinc_util.hh>

// External headers: BOINC
#ifdef BOINC
//#ifdef _WIN32
//#include <boinc_win.h>
//#endif
#include <boinc_api.h>

// C++ headers
#include <iostream>


namespace utility {
namespace boinc {

/// @brief Convert logical file names to physical names
void
resolve_filename( std::string & filename )
{
	using utility::file::file_exists;

	// first check if file exists
	if ( !file_exists( filename ) ) return;
	char resolved_name[ 512 ];
	if ( ::boinc_resolve_filename( filename.c_str(), resolved_name, sizeof( resolved_name ) ) ) {
		std::cout << "BOINC WARNING: can't resolve filename " << filename << std::endl;
	}
	filename = resolved_name;
}

} // namespace boinc
} // namespace utility

#endif
