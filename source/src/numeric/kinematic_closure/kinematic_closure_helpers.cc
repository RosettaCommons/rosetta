// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   kinematic_closure_helpers.cc
/// @brief  helpers for functions associated with bridgeObjects.cc
/// @author Daniel J. Mandell

// C++ headers
// Utility headers
#include <utility/vector1.hh>
// Rosetta Headers
#include <numeric/types.hh>
#include <iostream>
#include <iomanip>

#define SMALL 0.00001  // small epsilon to test for equality

using numeric::Real;

namespace numeric {
namespace kinematic_closure {

/// helper functions ///
void printVector(const utility::vector1<Real>& V) {
	for ( unsigned i=1; i<=V.size(); i++ ) {
		std::cout << V[i] << std::endl;
	}
	std::cout << std::endl;
}

/// @brief prints the matrix
/// @details This function used to intentionally print the transpose of the
/// matrix.  The rational was that "we use row-major indexing".  That didn't
/// make any sense to me, and I'd been mislead by the implicit transpose a
/// couple of times, so I got rid of it.
void printMatrix(const utility::vector1<utility::vector1<Real> >& M) {
	for ( unsigned i=1; i<=M.size(); i++ ) {
		for ( unsigned j=1; j<=M[1].size(); j++ ) {
			std::cout << std::setprecision(10) << std::setw(16) << M[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}

void printTranspose(const utility::vector1<utility::vector1<Real> >& M) {
	for ( unsigned i=1; i<=M[1].size(); i++ ) {
		for ( unsigned j=1; j<=M.size(); j++ ) {
			std::cout << std::setprecision(10) << std::setw(16) << M[j][i] << "\t";
		}
		std::cout << std::endl;
	}
}

// C is the product of matrices A and B. IE:
// C = A X B
void multMatrix(const utility::vector1<utility::vector1<Real> >& A,
	const utility::vector1<utility::vector1<Real> >& B,
	utility::vector1<utility::vector1<Real> >& C)
{
	unsigned cols, rows;
	cols = B.size();
	rows = B[1].size();
	C.resize(cols);
	for ( unsigned i=1; i<=cols; i++ ) {
		C[i].resize(rows);
		for ( unsigned j=1; j<=rows; j++ ) {
			C[i][j] = 0.0;
			for ( unsigned k=1; k<=rows; k++ ) {
				C[i][j] += A[k][j] * B[i][k];
			}
		}
	}
}

// C is the product of matrix A transposed and multiplied with matrix B. IE:
// C = A' X B
void multTransMatrix(const utility::vector1<utility::vector1<Real> >& A,
	const utility::vector1<utility::vector1<Real> >& B,
	utility::vector1<utility::vector1<Real> >& C)
{
	unsigned cols, rows;
	cols = B.size();
	rows = B[1].size();
	C.resize(cols);
	for ( unsigned i=1; i<=cols; i++ ) {
		C[i].resize(rows);
		for ( unsigned j=1; j<=rows; j++ ) {
			C[i][j] = 0.0;
			for ( unsigned k=1; k<=rows; k++ ) {
				C[i][j] += A[j][k] * B[i][k];
			}
		}
	}
}

// tests if all elements of two vectors are equal
bool vectorsEqual(const utility::vector1<Real>& A, const utility::vector1<Real>& B) {
	for ( unsigned i=1; i<=A.size(); i++ ) {
		if ( std::abs(long(A[i]) - long(B[i])) > SMALL ) {
			return false;
		}
	}
	return true;
}

} // end namespace kinematic_closure
} // end namespace numeric

