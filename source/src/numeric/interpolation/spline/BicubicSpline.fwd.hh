// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Forward declarations for the bicubic spline class.
/// @author Steven Combs
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_numeric_interpolation_spline_BicubicSpline_fwd_hh
#define INCLUDED_numeric_interpolation_spline_BicubicSpline_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace numeric {
namespace interpolation {
namespace spline {

class BicubicSpline;

typedef utility::pointer::shared_ptr< BicubicSpline > BicubicSplineOP;
typedef utility::pointer::shared_ptr< BicubicSpline const > BicubicSplineCOP;

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric

#endif /* BICUBIC_SPLINE_FWD_HH_ */
