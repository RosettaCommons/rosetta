// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathMatrix.srlz.hh
/// @brief  Serlialization routines for MathMatrix
/// @author Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_numeric_MathMatrix_srlz_HH
#define INCLUDED_numeric_MathMatrix_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <numeric/MathMatrix.hh>

namespace numeric {

template < class Archive, class T >
void
save_math_matrix( Archive & archive, numeric::MathMatrix< T > const & matrix )
{

	archive( matrix.get_number_rows() );
	archive( matrix.get_number_cols() );
	for ( platform::Size ii = 0; ii < matrix.get_number_rows(); ++ii ) {
		for ( platform::Size jj = 0; jj < matrix.get_number_cols(); ++jj ) {
			archive( matrix( ii, jj ) );
		}
	}

}

template < class Archive, class T >
void
load_math_matrix( Archive & archive, numeric::MathMatrix< T > & matrix )
{

	platform::Size nrows; archive( nrows );
	platform::Size ncols; archive( ncols );
	numeric::MathMatrix< T > second_matrix( nrows, ncols );
	for ( platform::Size ii = 0; ii < nrows; ++ii ) {
		for ( platform::Size jj = 0; jj < ncols; ++jj ) {
			archive( second_matrix(ii,jj) );
		}
	}
	matrix = second_matrix;
}

/// @brief Serialization function for arbitary MathMatrix classes.
template < class Archive, class T >
void save( Archive & archive, typename numeric::MathMatrix< T > const & matrix )
{
	save_math_matrix( archive, matrix );
}

/// @brief Deserialization function for MathMatrix's holding arbitrary data
template < class Archive, class T >
void load( Archive & archive, numeric::MathMatrix< T > & matrix )
{
	load_math_matrix( archive, matrix );
}


// template specialization for a handful of commonly used MathMatrixs; the implementation for these
// methods will live in MathMatrix.srlz.cc so that they don't have to be regenerated in as many
// object files.

/// @brief partial template specialization for MathMatrix of integers; the purpose of this specialization
/// is solely to reduce compile time and disk space, since MathMatrixs of integers are fairly common
template < class Archive > void save( Archive & archive, numeric::MathMatrix< float  > const & );
template < class Archive > void save( Archive & archive, numeric::MathMatrix< double > const & );

template < class Archive > void load( Archive & archive, numeric::MathMatrix< float  > & );
template < class Archive > void load( Archive & archive, numeric::MathMatrix< double > & );

}

#endif // SERIALIZATION

#endif // INCLUDED_numeric_MathMatrix_srlz_HH
