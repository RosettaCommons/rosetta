// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/antibody/grafting/RegExManager.fwd.hh
/// @brief   Forward declaration for RegExManager.
/// @author  Brian D. Weitzner <brian.weitzner@gmail.com>

#ifdef CXX11

#ifndef INCLUDED_protocols_antibody_grafting_regex_manager_FWD_HH
#define INCLUDED_protocols_antibody_grafting_regex_manager_FWD_HH

namespace protocols {
namespace antibody {
namespace grafting {

/// @brief  A singleton class for handling static const data for CDR-detecting regular expressions used by the RAb.
class RegExManager;

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif  // INCLUDED_protocols_antibody_grafting_regex_manager_FWD_HH

#endif // CXX11
