// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/linear_algebra/cholesky_decomposition.hh
/// @brief  Cholesky decomposition of a symmetric positive-definite matrix.
/// @author Andy Watkins

#ifndef INCLUDED_numeric_linear_algebra_cholesky_decomposition_hh
#define INCLUDED_numeric_linear_algebra_cholesky_decomposition_hh

#include <utility/vector1.hh>

namespace numeric {
namespace linear_algebra {

/// @brief Compute the lower-triangular Cholesky factor L of a symmetric positive-definite
/// covariance matrix, such that Sigma = L * L^T.
/// @details The input matrix must be square and symmetric positive-definite. Throws
/// utility::excn::Exception if the decomposition fails (e.g., the matrix is not
/// positive-definite, which can happen when user-specified correlation coefficients
/// are mutually inconsistent).
/// @param[in] covariance_matrix  Square symmetric positive-definite matrix, 1-indexed.
/// @return Lower-triangular matrix L (same dimensions as input), 1-indexed.
utility::vector1< utility::vector1< double > >
cholesky_factor(
	utility::vector1< utility::vector1< double > > const & covariance_matrix
);

} // namespace linear_algebra
} // namespace numeric

#endif // INCLUDED_numeric_linear_algebra_cholesky_decomposition_hh
