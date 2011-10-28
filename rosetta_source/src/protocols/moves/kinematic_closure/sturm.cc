// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   sturm.cc
/// @brief  implements sturm chain solver
/// @author Daniel J. Mandell

#include <protocols/moves/kinematic_closure/sturm.hh>


namespace protocols {
namespace moves {
namespace kinematic_closure {

double RELERROR;
int MAXIT, MAX_ITER_SECANT;

} // end namespace kinematic_closure
} // end namespace moves
} // end namespace protocols

