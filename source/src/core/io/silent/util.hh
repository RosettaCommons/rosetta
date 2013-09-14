// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/util.hh
///
/// @brief utility functions for silent-files
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_util_hh
#define INCLUDED_core_io_silent_util_hh

#include <map>
#include <string>

namespace core {
namespace io {
namespace silent {

/// @brief gzip all of the files in -out::file::silent().
void
gzip( void );

/////////////////////////////////////////////////////////////////
std::map< std::string, bool >
initialize_tag_is_done( std::string const & silent_file );

} // namespace silent
} // namespace io
} // namespace core

#endif
