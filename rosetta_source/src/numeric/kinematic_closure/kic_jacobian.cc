// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   kic_jacobian.cc
/// @brief  Jacobian for kinematic closure
/// @detail Adapted from Dodd, Boone & Theodorou, 1987 and Wu & Deem, 1999.                                              
/// @detail Need for random rotation is eliminated by combining rows 4 & 5                                              
/// @detail with correction factor to produce a 4x4 frame invariant form.                                                
/// @detail Uses formula eq.25, p.76, J.W. Gibbs, Vector Analysis, Dover, 1965  
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

#include <numeric/types.hh>
#include <utility/vector1.hh>

//#include <kinematic_closure_helpers.hh> // for printMatrix
#include <iostream>
#include <cmath>

using numeric::Real;
using numeric::Size;
using utility::vector1;
using std::cout;
using std::endl;

namespace numeric {
namespace kinematic_closure {

// prints the transpose of a matrix (since we store in column major format)
// delete and use kinematic_closure_helpers.hh version after integration
void printMatrix(const utility::vector1<utility::vector1<Real> >& M) {
	for (unsigned i=1; i<=M[1].size(); i++) {
		for (unsigned j=1; j<=M.size(); j++) {
			cout << M[j][i] << "\t";
		}
		cout << std::endl;
	}
	cout << endl;
}
	
// print a vector
// delete and use kinematic_closure_helpers.hh version after integration	
void printVector(const utility::vector1<Real>& V) {
	for (unsigned i=1; i<=V.size(); i++) {
		cout << V[i] << std::endl;
	}
	cout << std::endl;
}	
	
void printVectorTranspose(const utility::vector1<Real>& V) {
	for (unsigned i=1; i<=V.size(); i++) {
		cout << V[i] << "\t";
	}
	cout << std::endl;
}	
	
Real
det_3( const vector1< Real >& j1,
	   const vector1< Real >& j2,
	   const vector1< Real >& j3 )
{
	Real det_3, j11, j12, j13, j21, j22, j23, j31, j32, j33;
	
	j11 = j1[ 1 ];
	j12 = j1[ 2 ];
	j13 = j1[ 3 ];
	j21 = j2[ 1 ];
	j22 = j2[ 2 ];
	j23 = j2[ 3 ];
	j31 = j3[ 1 ]; 
	j32 = j3[ 2 ];
	j33 = j3[ 3 ];
	
	det_3 =   j11 * ( j22 * j33 - j23 * j32 )
	        - j12 * ( j21 * j33 - j23 * j31 )
	        + j13 * ( j21 * j32 - j22 * j31 );
	
	return det_3;
}

void
crossp( const vector1< Real >& U,
	    const vector1< Real >& V,
			  vector1< Real >& crossp_out )
{
	crossp_out[ 1 ] = U[ 2 ]*V[ 3 ] - U[ 3 ]*V[ 2 ];
	crossp_out[ 2 ] = U[ 3 ]*V[ 1 ] - U[ 1 ]*V[ 3 ];
	crossp_out[ 3 ] = U[ 1 ]*V[ 2 ] - U[ 2 ]*V[ 1 ];
}

Real
dot_product_3x3( const vector1< Real >& A,
				 const vector1< Real >& B ) {
	return A[ 1 ]*B[ 1 ] + A[ 2 ]*B[ 2 ] + A[ 3 ]*B[ 3 ];
}

void 
normalize( const vector1< Real >& V,
				 vector1< Real >& norm_out )
{
	Real norm = std::sqrt( dot_product_3x3( V, V ) );
	for( Size i=1; i<= 3; ++i ) {
		norm_out[ 1 ] = V[ 1 ] / norm;
		norm_out[ 2 ] = V[ 2 ] / norm;
		norm_out[ 3 ] = V[ 3 ] / norm;
	}
}

Real
kic_jac( const vector1< vector1< Real > >& r_n,
		 const vector1< vector1< Real > >& r_ca,
		 const vector1< vector1< Real > >& r_c )
{	
	vector1< vector1< Real > > axis( 6 ), R1( 6 ), R2( 6 ), j1_3( 4 );
	vector1< Real > R2_minus_R1( 3 ), norm_R2_minus_R1( 3 ), R1_6_minus_R1_i( 3 ), cross_axis_i_R1_6_minus_R1_i( 3 ),
		cross_axis_5_axis_6( 3 ), j4( 4 );
	Real det, jac;
	
	R1[ 1 ] = r_n [ 1 ]; R2[ 1 ] = r_ca [ 1 ];
	R1[ 2 ] = r_ca[ 1 ]; R2[ 2 ] = r_c  [ 1 ];
	R1[ 3 ] = r_n [ 2 ]; R2[ 3 ] = r_ca [ 2 ];
	R1[ 4 ] = r_ca[ 2 ]; R2[ 4 ] = r_c  [ 2 ];
	R1[ 5 ] = r_n [ 3 ]; R2[ 5 ] = r_ca [ 3 ];
	R1[ 6 ] = r_ca[ 3 ]; R2[ 6 ] = r_c  [ 3 ];
	
	//printMatrix( R1 );
	
	for( Size i=1; i<=6; ++i ) {
		for( Size j=1; j<=3; ++j ) {
			R2_minus_R1[ j ] = R2[ i ][ j ]-R1[ i ][ j ];
		}
		normalize( R2_minus_R1, norm_R2_minus_R1 );
		axis[ i ].resize( 3 );
		for ( Size j=1; j<=3; ++j ) {
			axis[ i ][ j ] = norm_R2_minus_R1[ j ];
		}
	}	
	//printMatrix( R1 );
	//printMatrix( axis );
	for( Size i=1; i<=4; ++i ) {
		j1_3[ i ].resize( 3 );
		// Jacobian elements (j11 - j14), (j21 - j24), (j31 - j34)
		for( Size j=1; j<=3; ++j ) {
			R1_6_minus_R1_i[ j ] = R1[ 6 ][ j ] - R1[ i ][ j ];
		}
		crossp( axis[ i ], R1_6_minus_R1_i, cross_axis_i_R1_6_minus_R1_i );
		//printVector( cross_axis_i_R1_6_minus_R1_i );
		for( Size j=1; j<=3; ++j ) {
			j1_3[ i ][ j ] = cross_axis_i_R1_6_minus_R1_i[ j ];
		}		
		// Jacobian elements (j41 - j45)
		crossp( axis[ 5 ], axis[ 6 ], cross_axis_5_axis_6 );
		j4[ i ] = dot_product_3x3( axis[ i ], cross_axis_5_axis_6 );
	}
	printMatrix( j1_3 );
	printVectorTranspose( j4 );
	
	det = -j4[ 1 ]*det_3( j1_3[ 2 ], j1_3[ 3 ], j1_3[ 4 ] )
		  +j4[ 2 ]*det_3( j1_3[ 1 ], j1_3[ 3 ], j1_3[ 4 ] )
		  -j4[ 3 ]*det_3( j1_3[ 1 ], j1_3[ 2 ], j1_3[ 4 ] )
		  +j4[ 4 ]*det_3( j1_3[ 1 ], j1_3[ 2 ], j1_3[ 3 ] );
	std::cout << det << std::endl;
	if( det == 0.0 ) {
		cout << "0.0 determinant; setting to 10^-100 to avoid zero division" << std::endl;
		det = std::pow( 10.0, -100.0);
		cout << "det: " << det << endl;
	}
	
	jac = 1.0 / std::abs( det );
	//cout << "Jacobian: " << jac << endl;
	return jac;
}

void
test_kic_jac() {
	Real r_n_vals[]  = { 3.289, 1.561, 0.295, 15.226, -2.198, -11.850, 15.763, 1.816, 2.790 };
	Real r_ca_vals[] = { 3.932, 2.871, 0.330, 16.145, -1.322, -12.567, 15.339, 2.953, 1.982 };
	Real r_c_vals[]  = { 4.604, 3.183, -1.003, 17.593, -1.638, -12.201, 13.822, 2.976, 1.830 };
	Real jac;
	vector1< vector1< Real > > r_n( 3 ), r_ca( 3 ), r_c( 3 );
	Size N=3;
		
	for( Size i=1; i<=3; ++i ) {
		r_n [ i ].resize( 3 );
		r_ca[ i ].resize( 3 );
		r_c [ i ].resize( 3 );
	}
		
	Size ind=0;
	for( Size i=1; i<=3; ++i ) {
		for( Size j=1; j<=N; ++j ) {
			r_n [ i ][ j ] = r_n_vals [ ind ];
			r_ca[ i ][ j ] = r_ca_vals[ ind ];
			r_c [ i ][ j ] = r_c_vals [ ind ];
			++ind;
		}
	}
		
	printMatrix( r_n );
	//printMatrix( r_ca );
	//printMatrix( r_c );
	
	jac = kic_jac( r_n, r_ca, r_c );
	
	cout << "kic_jac: " << jac << endl;
	// output should be "kic_jac: 0.000623412"
}

void
test_normalize() {
	vector1< Real > V( 3 ), V_norm( 3 );
	
	V[ 1 ] = 0.64299988746643066;
	V[ 2 ] = 1.3100000619888306;
	V[ 3 ] = 0.03500002622604370;
	
	normalize( V, V_norm );
	
	cout << "V_norm: ";
	for( Size i=1; i<=3; ++i ) {
		cout << V_norm[ i ] << "\t";
	}
	cout << endl;
	
	// output should be "V_norm: 0.440496        0.897434        0.0239773"
}

void
test_crossp() {
	vector1< Real > U( 3 ), V( 3 ), crossp_out( 3 );
	
	U[ 1 ] = 0.44049623473987115;
	U[ 2 ] = 0.89743420809728203;
	U[ 3 ] = 0.02397726666658848;
	V[ 1 ] = 12.049999713897705;
	V[ 2 ] = 1.3920000791549683;
	V[ 3 ] = 1.6870000064373016;
	
	crossp( U, V, crossp_out );
	
	cout << "crossp_out: ";
	for( Size i=1; i<=3; ++i ) {
		cout << crossp_out[ i ] << "\t";
	}
	cout << endl;
	// output should be "crossp_out: 1.4806      -0.454191       -10.2009"
}

void
test_det_3() {	
	vector1< Real > j1( 3 ), j2( 3 ), j3( 3 );
	Real det_3_out;
	
	j1[ 1 ] = 0.40964212899059071;
	j1[ 2 ] = -10.698361160396704;
	j1[ 3 ] = -2.2975314412577430;
	j2[ 1 ] = 10.843018060354922;
	j2[ 2 ] = -8.7735573087774217;
	j2[ 3 ] = 3.1786677686625100;
	j3[ 1 ] = -4.0365017319290022;
	j3[ 2 ] = -13.993128305123419;
	j3[ 3 ] = 3.8880472947998008;	
		
	cout << "j1:\t";
	for ( Size i = 1; i <= 3; ++i ) {
		cout << j1[ i ] << "\t";
	}
	cout << endl;
	
	cout << "j2:\t";
	for ( Size i = 1; i <= 3; ++i ) {
		cout << j2[ i ] << "\t";
	}
	cout << endl;
	
	cout << "j3:\t";
	for ( Size i = 1; i <= 3; ++i ) {
		cout << j3[ i ] << "\t";
	}
	cout << endl;	
	
	det_3_out = det_3( j1, j2, j3 );
	cout << "det_3: " << det_3_out << endl;
	// output should be "det_3: 1022.5"
}
	
} // end namespace kinematic_closure
} // end namespace numeric

int
main( int argc, char ** argv )
{
	//numeric::kinematic_closure::test_det_3();
	//numeric::kinematic_closure::test_crossp();
	//numeric::kinematic_closure::test_normalize();
	numeric::kinematic_closure::test_kic_jac();
	
}
