// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/boinc/boinc_util.hh
/// @brief  Wrappers to make BOINC work
/// @author David Kim (dekim@u.washington.edu)


#ifndef INCLUDED_utility_boinc_boinc_util_hh
#define INCLUDED_utility_boinc_boinc_util_hh

#ifdef BOINC

// C++ headers
#include <string>

namespace utility {
namespace boinc {

/// @brief Convert logical file names to physical names
void
resolve_filename( std::string & filename );

} // namespace boinc
} // namespace utility

#endif


#endif // INCLUDED_utility_boinc_boinc_util_HH
