// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/Interpolator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///

#ifndef INCLUDED_numeric_interpolation_spline_Interpolator_hh
#define INCLUDED_numeric_interpolation_spline_Interpolator_hh

#include <numeric/types.hh>

#include <utility/vector1.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

namespace numeric {
namespace interpolation {
namespace spline {

using numeric::Real;

class Interpolator : public utility::pointer::ReferenceCount {

public:

	virtual void interpolate( Real x, Real & y, Real & dy ) = 0;

};

typedef utility::pointer::owning_ptr< Interpolator > InterpolatorOP;

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
