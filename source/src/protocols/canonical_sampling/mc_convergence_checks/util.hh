// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_util_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_util_hh


// type headers

// unit headers
#include <protocols/moves/MonteCarlo.fwd.hh>

// package headers

// utility headers
// #include "utility/basic_sys_util.h"

// C++ headers

// Forward declarations

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

extern void setup_convergence_checks_from_cmdline( moves::MonteCarlo& mc );

} // mc_convergence_check
} // moves
} // rosetta

#endif
