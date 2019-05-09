// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/linear_algebra/minimum_bounding_ellipse.cc
/// @brief Given a set of points, calculate the minimum bounding ellipse
/// @author Rebecca Alford (ralford3@jhu.edu)

// Unit Headers
#include <numeric/linear_algebra/minimum_bounding_ellipse.hh>

// Project Headers
#include <numeric/MathMatrix.hh>
#include <numeric/MathMatrix_operations.hh>
#include <numeric/xyzVector.hh>
#include <numeric/linear_algebra/singular_value_decomposition.hh>
#include <numeric/linear_algebra/EllipseParameters.hh>

#include <numeric/types.hh>

// Basic/utility Headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C++ Headers
#include <cstdlib>
#include <cmath>

namespace numeric {
namespace linear_algebra {


/// @brief Use the Khachiyan Algorithm to compute the minimum volume enclosing ellipsoid given a set of (x,y) data points
/// @note Implementation taken from Stack Overflow, adapted from MATLAB minimum volume enclosing ellipsoid code by Nima Moshtagh
/// @note https://stackoverflow.com/questions/1768197/bounding-ellipse
EllipseParametersOP
minimum_bounding_ellipse(
	utility::vector1< xyzVector< Real > > points,
	Real tolerance,
	Size max_iterations
) {

	using namespace numeric;
	using namespace utility;

	// Initialize constant dimensions
	Size const dimension( 2 );
	Size const N( points.size() );

	// create matrix P from the points vector
	MathMatrix< Real > P = MathMatrix< Real >( 2, points.size() );
	for ( Size ii = 0; ii < points.size(); ++ii ) {
		P( 0, ii ) = points[ii+1].x();
		P( 1, ii ) = points[ii+1].y();
	}

	// create matrix Q from matrix P
	MathVector< Real > Q_col2 = MathVector< Real >( 1, points.size() );
	MathMatrix< Real > Q = MathMatrix< Real >( 3, points.size() );
	for ( Size ii = 0; ii < points.size(); ++ii ) {
		Q( 0, ii ) = points[ii+1].x();
		Q( 1, ii ) = points[ii+1].y();
		Q( 2, ii ) = 1;
	}

	// Initialize counts and errors
	Size count( 1 );
	Real error( 1 );

	// Initialize u vector: An Nx1 vector where each element is 1/N
	Real const fill_value = 1/static_cast< Real >(N);
	MathMatrix< Real > u = MathMatrix< Real >( N, 1, fill_value );

	// One initial difference between the C++ and MATLAB code is that
	// Khachiyan Algorithm
	while ( error > tolerance ) {

		// place the elements of u in an NxN vector
		MathMatrix< Real > diag_u = MathMatrix< Real >( N, N );
		for ( Size ii = 0; ii <= N-1; ++ii ) {
			diag_u( ii, ii ) = u(ii, 0);
		}

		// multiply Q * diag(u) * Q' (Q' is the transpose)
		MathMatrix< Real > X = (Q * diag_u) * non_square_transpose(Q);
		MathMatrix< Real > X_inverse = X.inverse_square_matrix();
		MathMatrix< Real > M_prime = (non_square_transpose(Q) * X_inverse) * Q;

		// Get the diagonal of M prime
		utility::vector1< Real > M = utility::vector1< Real>(N);
		for ( Size ii = 0; ii <= N-1; ++ii ) {
			M[ii+1] = M_prime(ii, ii);
		}

		// Find the value and location of the maximum element
		Real maximum = max(M);
		Size j = arg_max(M);

		// Calculate the ascent step size
		Real step_size = (maximum - dimension - 1)/((dimension+1)*(maximum-1));

		// Calculate a new u vector
		MathMatrix< Real > new_u = MathMatrix< Real >(N, 1);
		for ( Size ii = 0; ii < N; ++ii ) {
			new_u(ii,0) = u(ii,0) * (1-step_size);
		}
		new_u(j-1,0) = new_u(j-1,0) + step_size;

		// Compute error as the square root of the SSD between new_u and uf
		error = sum_of_square_differences( u, new_u );

		// Increment count and replace u
		count = count + 1;
		u = new_u;

		// If at max iterations - go one more iteration then stop
		if ( count == max_iterations ) {
			error = tolerance;
		}
	}

	// Place elements in the u vector into the diagonal of an NxN matrix
	numeric::MathMatrix< Real > U = numeric::MathMatrix< Real >(N, N);
	for ( Size ii = 0; ii < N; ++ii ) {
		U(ii, ii) = u(ii,0);
	}

	// Compute the A-matrix with parameters
	MathMatrix< Real > c = P * u;
	MathMatrix< Real > A_prime1 = (P * U) * non_square_transpose(P);
	MathMatrix< Real > A_prime2 = c*non_square_transpose(c);
	Real inv = 1/static_cast< Real >(dimension);
	MathMatrix< Real > A = (A_prime1 - A_prime2).inverse() * inv;

	// Convert the A matrix to vector1< vector1< Real > > format
	vector1< vector1< Real > > A_converted;
	A_converted.resize( A.get_number_rows() );
	for ( Size ii = 0; ii < A.get_number_rows(); ++ii ) {
		for ( Size jj = 0; jj < A.get_number_cols(); ++jj ) {
			A_converted[ii+1].push_back( A(ii, jj) );
		}
	}

	// The A matrix should be dimension x dimension (2x2 in this case
	vector1< Real > W( dimension, 0 );
	vector1< vector1< Real > > V( dimension, vector1< Real >( dimension, 0 ) ); // n x n matrix
	linear_algebra::svdcmp( A_converted, dimension, dimension, W, V );

	// Convert the rotation matrix V to a MathMatrix
	MathMatrix< Real > rotation = MathMatrix< Real >(2,2);
	for ( Real ii = 1; ii <= V.size(); ++ ii ) {
		for ( Real jj = 1; jj <= V[ii].size(); ++jj ) {
			rotation(ii-1,jj-1) = V[ii][jj];
		}
	}

	// Calculate the center from the C matrix
	Real h( c(0,0) );
	Real k( c(1,0) );
	Real a = 1/sqrt( W[1] );
	Real b = 1/sqrt( W[2] );

	EllipseParametersOP ellipse( new EllipseParameters( h, k, a, b, rotation ) );
	return ellipse;

}


/// @brief Calculate the sum-of-square differences between
/// values stored in two vector1 objects
Real
sum_of_square_differences(
	MathMatrix< Real > old_u,
	MathMatrix< Real > new_u
) {

	// Calculate the difference between the two matrices
	MathMatrix< Real > diff = new_u - old_u;
	Real a(0);
	for ( Size ii = 0; ii < diff.get_number_rows(); ++ii ) {
		a = a + pow( std::abs(diff(ii,0)), 2 );
	}
	return sqrt(a);
}

/// @brief Calculate the transpose of a non-square MathMatrix
/// and return the result as a new MathMatrix
MathMatrix< Real >
non_square_transpose(
	MathMatrix< Real > matrix_in
) {

	MathMatrix< Real > matrix_out = MathMatrix< Real >( matrix_in.get_number_cols(), matrix_in.get_number_rows() );

	for ( Size ii = 0; ii < matrix_in.get_number_cols(); ++ii ) {
		for ( Size jj = 0; jj < matrix_in.get_number_rows(); ++jj ) {
			matrix_out(ii,jj) = matrix_in(jj,ii);
		}
	}

	return matrix_out;
}

/// @brief Check whether a given test point lies within an ellipse
bool
point_in_ellipse(
	xyzVector< Real > p,
	Real const h,
	Real const k,
	Real const a,
	Real const b,
	MathMatrix< Real > rotation
) {

	Real first_term = (rotation(0,0)*(p.x()-h) + rotation(1,0)*(p.y()-k))/a;
	Real second_term = (rotation(1,0)*(p.x()-h) - rotation(0,0)*(p.y()-k))/b;
	Real lhs = pow( first_term, 2 ) + pow( second_term, 2 );

	if ( lhs <= 1 ) {
		return true;
	} else {
		return false;
	}

}


} // ns linear_algebra
} // ns numeric

