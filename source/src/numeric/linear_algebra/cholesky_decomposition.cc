// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/linear_algebra/cholesky_decomposition.cc
/// @brief  Cholesky decomposition of a symmetric positive-definite matrix.
/// @author Andy Watkins

#include <numeric/linear_algebra/cholesky_decomposition.hh>

#include <utility/excn/Exceptions.hh>

#include <Eigen/Cholesky>

namespace numeric {
namespace linear_algebra {

utility::vector1< utility::vector1< double > >
cholesky_factor(
	utility::vector1< utility::vector1< double > > const & covariance_matrix
) {
	auto const n = covariance_matrix.size();
	if ( n == 0 ) {
		throw CREATE_EXCEPTION( utility::excn::Exception,
			"Error in numeric::linear_algebra::cholesky_factor(): Empty matrix provided." );
	}

	for ( platform::Size i = 1; i <= n; ++i ) {
		if ( covariance_matrix[i].size() != n ) {
			throw CREATE_EXCEPTION( utility::excn::Exception,
				"Error in numeric::linear_algebra::cholesky_factor(): Matrix is not square." );
		}
	}

	Eigen::MatrixXd mat( n, n );
	for ( platform::Size i = 1; i <= n; ++i ) {
		for ( platform::Size j = 1; j <= n; ++j ) {
			mat( i - 1, j - 1 ) = covariance_matrix[i][j];
		}
	}

	Eigen::LLT< Eigen::MatrixXd > llt( mat );
	if ( llt.info() != Eigen::Success ) {
		throw CREATE_EXCEPTION( utility::excn::Exception,
			"Error in numeric::linear_algebra::cholesky_factor(): Cholesky decomposition failed. "
			"The covariance matrix is not positive-definite. This typically means the specified "
			"correlation coefficients are mutually inconsistent (e.g., A correlates strongly with "
			"B and C, but B and C anti-correlate in a way that is mathematically impossible)." );
	}

	Eigen::MatrixXd const L = llt.matrixL();

	utility::vector1< utility::vector1< double > > result( n, utility::vector1< double >( n, 0.0 ) );
	for ( platform::Size i = 1; i <= n; ++i ) {
		for ( platform::Size j = 1; j <= n; ++j ) {
			result[i][j] = L( i - 1, j - 1 );
		}
	}

	return result;
}

} // namespace linear_algebra
} // namespace numeric
