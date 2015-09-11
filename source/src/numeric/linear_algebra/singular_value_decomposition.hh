// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/linear_algebra/singular_value_decomposition.hh
/// @brief  Header file for Singular Value Decomposition.

#ifndef INCLUDED_numeric_linear_algebra_singular_value_decomposition_HH
#define INCLUDED_numeric_linear_algebra_singular_value_decomposition_HH

#include <numeric/types.hh>
#include <utility/vector1.fwd.hh>

namespace numeric {
namespace linear_algebra {

void
svdcmp( utility::vector1< utility::vector1< Real > > &a,
	    Size const m, Size const n, 
	    utility::vector1< Real > &w,
	    utility::vector1< utility::vector1< Real > > &v );


} // end namespace linear_algebra
} // end namespace numeric

#endif // INCLUDED_numeric_linear_algebra_singular_value_decomposition_HH
