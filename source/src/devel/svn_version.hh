// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/svn_version.hh
///
/// @brief
/// @author Ian W. Davis
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_devel_svn_version_hh
#define INCLUDED_devel_svn_version_hh

#include <string>

namespace devel {

// @brief This function, called in devel::init(), initializes two peices of data in the core library
// indicating the version of the code based on the SVN or git repository the executable was compiled from
void register_version_with_core();

} // namespace devel

#endif // INCLUDED_devel_svn_version_HH
