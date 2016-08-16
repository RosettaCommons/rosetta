// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/cyclic_coordinate_descent.hh
/// @brief Compute the angle that minimizes the deviation between a set of atoms
/// @author Brian D. Weitzner
/// @author Labonte <JWLabonte@jhu.edu>
/// @author Phil Bradley

#ifndef INCLUDED_numeric_cyclic_coordinate_descent_HH
#define INCLUDED_numeric_cyclic_coordinate_descent_HH

#include <numeric/types.hh>

#include <utility/vector1.hh>

namespace numeric {

void ccd_angle(
	utility::vector1< xyzVector< Real > > const & F,
	utility::vector1< xyzVector< Real > > const & M,
	xyzVector< Real > const & axis_atom,
	xyzVector< Real > const & theta_hat,
	Real & alpha,
	Real & S );

} // numeric

#endif // INCLUDED_numeric_cyclic_coordinate_descent_HH
