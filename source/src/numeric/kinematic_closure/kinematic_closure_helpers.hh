// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   kinematic_closure_helpers.hh
/// @brief  Header file for kinematic_closure_helpers.cc
/// @author Daniel J. Mandell

#ifndef INCLUDED_numeric_kinematic_closure_kinematic_closure_helpers_hh
#define INCLUDED_numeric_kinematic_closure_kinematic_closure_helpers_hh

// Rosetta Headers
#include <numeric/types.hh>

// Utility headers

#include <utility/vector1.fwd.hh>
#include <utility/fixedsizearray1.fwd.hh>

#include <iomanip>

namespace numeric {
namespace kinematic_closure {

/// helper functions ///
template< class Iterable >
void printVector( const Iterable & V ) {
	for ( auto const elem : V ) {
		std::cout << elem << std::endl;
	}
	std::cout << std::endl;
}
	
/// @brief prints the matrix
/// @details This function used to intentionally print the transpose of the
/// matrix.  The rational was that "we use row-major indexing".  That didn't
/// make any sense to me, and I'd been mislead by the implicit transpose a
/// couple of times, so I got rid of it.
template< class MatrixLike >
void printMatrix(const MatrixLike & M) {
	for ( auto const & outer : M ) {
		for ( auto const & elem : outer ) {
			std::cout << std::setprecision(10) << std::setw(16) << elem << "\t";
		}
		std::cout << std::endl;
	}
}
	
template< class MatrixLike >
void printTranspose( const MatrixLike & M) {
	for ( unsigned i=1; i<=M[1].size(); i++ ) {
		for ( unsigned j=1; j<=M.size(); j++ ) {
			std::cout << std::setprecision(10) << std::setw(16) << M[j][i] << "\t";
		}
		std::cout << std::endl;
	}
}

template< Size rows >
void multMatrix(const utility::fixedsizearray1<utility::fixedsizearray1<numeric::Real, rows >, rows >& A,
	const utility::vector1<utility::fixedsizearray1<numeric::Real, rows> >& B,
				utility::vector1<utility::fixedsizearray1<numeric::Real, rows> >& C)
{
	Size const cols = B.size();
	C.resize( cols );
	for ( unsigned i=1; i<=cols; i++ ) {
		for ( unsigned j=1; j<=rows; j++ ) {
			C[i][j] = 0.0;
			for ( unsigned k=1; k<=rows; k++ ) {
				C[i][j] += A[k][j] * B[i][k];
			}
		}
	}
}

template< Size rows >
void multTransMatrix(const utility::fixedsizearray1<utility::fixedsizearray1<numeric::Real, rows>, rows >& A,
	const utility::vector1<utility::fixedsizearray1<numeric::Real, rows> >& B,
					 utility::vector1<utility::fixedsizearray1<numeric::Real, rows> >& C) {
	Size const cols = B.size();
	C.resize( cols );
	for ( unsigned i=1; i<=cols; i++ ) {
		for ( unsigned j=1; j<=rows; j++ ) {
			C[i][j] = 0.0;
			for ( unsigned k=1; k<=rows; k++ ) {
				C[i][j] += A[j][k] * B[i][k];
			}
		}
	}
}
	
bool vectorsEqual(const utility::vector1<numeric::Real>& A, const utility::vector1<numeric::Real>& B);

// Do not call. By the way, there's a better way to do this but I forget...
void force_instantiation();

} // end namespace kinematic_closure
} // end namespace numeric

#endif
