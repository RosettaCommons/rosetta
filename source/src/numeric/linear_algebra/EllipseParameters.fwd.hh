// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/linear_algebra/EllipseParameters.fwd.hh
/// @brief Container class for ellipse parameters
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_numeric_linear_algebra_EllipseParameters_fwd_hh
#define INCLUDED_numeric_linear_algebra_EllipseParameters_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace numeric {
namespace linear_algebra {

class EllipseParameters;

typedef utility::pointer::shared_ptr< EllipseParameters > EllipseParametersOP;
typedef utility::pointer::shared_ptr< EllipseParameters const > EllipseParametersCOP;

} //numeric
} //linear_algebra

#endif //INCLUDED_numeric_linear_algebra_EllipseParameters_fwd_hh
