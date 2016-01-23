// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/util.cc
/// @brief Helpers for antibody grafting code
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__
#include <regex>
#include <string>
#endif

namespace protocols {
namespace antibody {
namespace grafting {

// It's only usable if the regex works.
bool antibody_grafting_usable() {
#ifndef __ANTIBODY_GRAFTING__
	return false;
#else
	std::string s("is the regex code working?");
	std::regex reg("regex");
	return std::regex_search(s,reg); // Problematic GCC versions return false for everything
#endif
}

} // grafting
} // antibody
} // protocols
