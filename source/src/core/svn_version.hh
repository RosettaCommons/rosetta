// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/svn_version.hh
///
/// @brief
/// @author Ian W. Davis
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_svn_version_hh
#define INCLUDED_core_svn_version_hh

#include <string>

namespace core {

/// @brief Initialize this data from a lower-level library at startup to avoid relinking every time
/// you update your svn / git version.  This function should be called at most once.
void set_svn_version_and_url( std::string const & version, std::string const & url );

/// @brief Read access to the svn / git version.
std::string minirosetta_svn_version();

/// @brief Read access to the svn / git url.
std::string minirosetta_svn_url();

} // namespace core

#endif // INCLUDED_core_svn_version_HH
