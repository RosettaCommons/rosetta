// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// TODO: Get rid of this class

/// @file  core/conformation/membrane/util.fwd.hh
/// @brief  utility functions for membrane things
/// @author  JKLeman (julia.koehler1982@gmail.com)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_util_hh
#define INCLUDED_core_conformation_membrane_util_hh

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>

// Package Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

namespace core {
namespace conformation {
namespace membrane {

// read multiple spanfiles from options
utility::vector1< std::string > spanfile_names();

// read single spanfile from options
std::string spanfile_name();

}// membrane
}// conformation
}// core

#endif // INCLUDED_core_conformation_membrane_util_hh
