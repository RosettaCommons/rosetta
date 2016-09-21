// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/interpolation/spline/PolycubicSplineBase.fwd.hh.
/// @brief Forward declaration of PolycubicSplineBase class, and declaration of (templated) owning pointers.
/// @details Note that, if you declare an owning pointer to a PolycubicSpline, you must know the dimensionality
/// at pointer declaration time!  If you want a general owning pointer to an arbitrary PolycubicSpline of
/// dimensionality to be determined later, use an owning pointer to the PolycubicSplineBase class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_numeric_interpolation_spline_PolycubicSplineBase_fwd_hh
#define INCLUDED_numeric_interpolation_spline_PolycubicSplineBase_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace numeric {
namespace interpolation {
namespace spline {

// Forward declaration:
class PolycubicSplineBase;

// Owning pointer:
typedef utility::pointer::shared_ptr< PolycubicSplineBase > PolycubicSplineBaseOP;

// Const-access owning pointer:
typedef utility::pointer::shared_ptr< PolycubicSplineBase const > PolycubicSplineBaseCOP;

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric


#endif /* INCLUDED_numeric_interpolation_spline_PolycubicSplineBase_fwd_hh */
