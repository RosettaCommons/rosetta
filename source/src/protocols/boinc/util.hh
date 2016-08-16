// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/boinc/util
/// @brief utility functions for filtering boinc outputs
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_boinc_util_hh
#define INCLUDED_protocols_boinc_util_hh

#include <map>
#include <string>
#include <core/types.hh>

namespace protocols {
namespace boinc {

/////////////////////////////////////////////////////////////////
void boincOutputFilter(core::Real runTime, core::Real minTimePerModel);

} // namespace protocols
} // namespace boinc

#endif
