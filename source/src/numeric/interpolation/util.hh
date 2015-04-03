// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/util.hh
/// @author Sam DeLuca
/// @author Stephanie Hirst

#ifndef INCLUDED_numeric_interpolation_util_HH
#define INCLUDED_numeric_interpolation_util_HH

#include <numeric/interpolation/spline/SplineGenerator.hh>

namespace numeric {
namespace interpolation {

/// @brief given a file, return a 2D spline
spline::SplineGenerator spline_from_file(std::string const & filename, platform::Real const & bin_size);
}
}

#endif /* INCLUDED_numeric_interpolation_util_HH */
