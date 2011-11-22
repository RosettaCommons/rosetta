// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNAResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/determinants.hh>
//////////////////////////////////
#include <core/types.hh>

#include <iostream>
// AUTO-REMOVED #include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <set>

//Auto Headers
#include <utility/vector1.hh>
using namespace core;

namespace protocols {
namespace swa {
namespace rna {

///////////////////////////////////////////////////////////////////////////////
void
output_matrix( utility::vector1< utility::vector1< Real > > const M ){
	for( Size i = 1; i <= M.size(); i++ ) {
		for( Size j = 1; j <= M[i].size(); j++ ) {
			std::cout << ' ' << ObjexxFCL::fmt::F(8,4,M[i][j]);
		}
		std::cout << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
Real
get_determinant( utility::vector1< utility::vector1< Real > > const M ){

	utility::vector1< utility::vector1< Real > > A = M;

	// Nice, geometric algorithm from detm.m (matlab )
	Size const m = A.size();
	Size const n = A[1].size();

	Real d( 1.0 );

	for (Size k = 1; k <= n; k++ ) {

		// Is this swapping step necessary?
		Size max_element( 0 );
		Real max_val( 0.0 );
		for ( Size h = k; h <= m; h++ ){
			if ( std::abs( A[h][k] ) > max_val ){
				max_val = std::abs( A[h][k] );
				max_element = h;
			}
		}
		if ( max_element == 0 ) return 0.0;
		Size const h = max_element;

		if ( k != h ){
			//swap
			utility::vector1< Real > x = A[k];
			A[k] = A[h];
			A[h] = x;
			d *= -1.0;
		}

		Real const c = A[k][k];
		//		A[ k ] = A[ k ]/c; //scale
		for ( Size j = 1; j <= n; j++ ) A[k][j] = A[k][j] / c;
		d = c * d;

		for ( Size i = (k+1); i <= m; i++ ){
			// A[i] = A[i] - A[i][k]*A[k]; // shear
			Real const save_coefficient = A[i][k];
			for ( Size j = 1; j <= n; j++ ) {
				A[i][j] = A[i][j] - save_coefficient * A[k][j];
			}
		}

		//std::cout << " AFTER " << k << std::endl;
		//output_matrix( A );

	}

	return d;

}

///////////////////////////////////////////////////////////////
// Real
// get_determinant_armadillo( utility::vector1< utility::vector1< Real > > const & J ){
// 	using namespace arma;
// 	mat m;
// 	m.zeros( 6, 6 );

// 	//	std::cout << "TAKING DETERMINANT OF MATRIX: " << std::endl;
// 	for ( Size i = 0; i < 6; i++ ) {
// 		for ( Size j = 0; j < 6; j++ ) {
// 			m( i, j ) = J[i+1][j+1];
// 			//			std::cout << ' ' << m( i, j );
// 		}
// 		//		std::cout << std::endl;
// 	}

// 	Real const determinant = det( m );
// 	//	std::cout << "   determinant: " << determinant << std::endl;
// 	return determinant;
// }





}
}
}
