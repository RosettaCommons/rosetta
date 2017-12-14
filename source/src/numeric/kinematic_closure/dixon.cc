// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   dixon.cc
/// @brief  computes the dixon resultant
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

// C++ headers
// Utility headers
//#include <basic/Tracer.hh> // tracer output

// Unit headers
#include <numeric/kinematic_closure/dixon.hh>

// Numeric headers
#include <numeric/types.hh>
#include <numeric/kinematic_closure/sturm.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/linear_algebra/GeneralizedEigenSolver.hh>
#include <utility/fixedsizearray1.hh>

// External headers
#include <Eigen/Dense>

// C++ headers
#include <cmath>
#include <iostream>

// Constants
#define PP4x2_VECSIZE 7
#define DIXON_SIZE 8
#define DIXON_RESULTANT_SIZE 3

using numeric::Real;
using utility::vector1;

namespace numeric {
namespace kinematic_closure {

using namespace std;
Eigen::IOFormat scipy(8, 0, ", ", ",\n", "[", "]", "[", "]");
string skipl = "\n\n";

typedef Eigen::Matrix<Real, 8, 8> Matrix8;
typedef Eigen::Matrix<Real, 16, 16> Matrix16;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1, 0, 16, 1> Vector16;
using SolverType = linear_algebra::GeneralizedEigenSolver<Matrix16>;

// If I can rant for a little bit, using STL vectors as a matrix type was a
// fucking horrible idea.  For one thing, the API presented by the std::vector
// class is totally inappropriate for doing linear algebra.  More importantly,
// you run the risk of getting really terrible caching performance because
// these matrices aren't necessarily stored contiguously in memory.

// AMW: changing to fixedsizearray1
typedef utility::fixedsizearray1<Real,3> PseudoVector;
typedef utility::fixedsizearray1<PseudoVector,3> PseudoMatrix;

/// dixon functions ///

// point_value2 {{{1
/* C becomes the point value 2 between A and t */
void point_value2(const utility::fixedsizearray1<Real,3>& A, const utility::vector1<Real>& t, utility::vector1<Real>& C) {
	C.resize(t.size());
	for ( unsigned i=1; i<=t.size(); i++ ) {
		C[i] = A[1]+t[i] * (A[2]+t[i]*A[3]);
	}
}

// point_value4 {{{1
/* C becomes the point value 4 between A and t */
void point_value4(const utility::vector1<Real>& A, const utility::vector1<Real>& t, utility::vector1<Real>& C) {
	C.resize(t.size());
	for ( unsigned i=1; i<=t.size(); i++ ) {
		C[i] = A[1]+t[i] * (A[2] + t[i] * (A[3] + t[i] * (A[4] + t[i]*A[5])));
	}
}

// point_value6 {{{1
/* C becomes the point value 6 between B and t */
void point_value6(const utility::vector1<Real>& B, const utility::vector1<Real>& t, utility::vector1<Real>& C) {
	C.resize(t.size());
	for ( unsigned i=1; i<=t.size(); i++ ) {
		C[i] = B[1]+t[i]*(B[2]+t[i]*(B[3]+t[i]*(B[4]+t[i]*(B[5]+t[i]*(B[6]+t[i]*B[7])))));
	}
}

// point_value8 {{{1
/* C becomes the point value 8 between A and t */
void point_value8(const utility::vector1<Real>& A, const utility::vector1<Real>& t, utility::vector1<Real>& C) {
	C.resize(t.size());
	for ( unsigned i=1; i<=t.size(); i++ ) {
		C[i] = A[1]+t[i]*(A[2]+t[i]*(A[3]+t[i]*(A[4]+t[i]*(A[5]+t[i]*(A[6]+t[i]*(A[7]+t[i]*(A[8]+t[i]*A[9])))))));
	}
}

// point_value16 {{{1
/* C becomes the point value 16 between A and t */
void point_value16(const utility::vector1<Real>& A, const utility::vector1<Real>& t, utility::vector1<Real>& C) {
	C.resize(t.size());
	for ( unsigned i=1; i<=t.size(); i++ ) {
		C[i] = A[1]+t[i]*(A[2]+t[i]*(A[3]+t[i]*(A[4]+t[i]*(A[5]+t[i]*(A[6]+t[i]*(A[7]+t[i]*(A[8]+t[i]*(A[9]+t[i]*(A[10]+t[i]*(A[11]+t[i]*(A[12]+t[i]*(A[13]+t[i]*(A[14]+t[i]*(A[15]+t[i]*(A[16]+t[i]*A[17])))))))))))))));
	}
}
// }}}1

// polyProduct2x2 {{{1
/* C becomes the polyProduct2x2 of A and B */
void polyProduct2x2(const utility::fixedsizearray1<Real,3>& A, const utility::fixedsizearray1<Real,3>& B, utility::vector1<Real>& C) {
	C.resize(PP4x2_VECSIZE);
	C[5] = A[3]*B[3];
	C[4] = A[3]*B[2]+A[2]*B[3];
	C[3] = A[3]*B[1]+A[2]*B[2]+A[1]*B[3];
	C[2] = A[2]*B[1]+A[1]*B[2];
	C[1] = A[1]*B[1];
	return;
}

// polyProduct4x2 {{{1
/* C becomes the polyProduct4x2 of A and B */
void polyProduct4x2(const utility::vector1<Real>& A, const utility::fixedsizearray1<Real,3>& B, utility::vector1<Real>& C) {
	C.resize(7);
	C[7] = A[5]*B[3];
	C[6] = A[5]*B[2]+A[4]*B[3];
	C[5] = A[5]*B[1]+A[4]*B[2]+A[3]*B[3];
	C[4] = A[4]*B[1]+A[3]*B[2]+A[2]*B[3];
	C[3] = A[3]*B[1]+A[2]*B[2]+A[1]*B[3];
	C[2] = A[2]*B[1]+A[1]*B[2];
	C[1] = A[1]*B[1];
	return;
}

// polyProduct4x4 {{{1
/* C becomes the polyProduct4x4 of A and B */
void polyProduct4x4(const utility::vector1<Real>& A, const utility::vector1<Real>& B, utility::vector1<Real>& C) {
	C.resize(9);
	C[9] = A[5]*B[5];
	C[8] = A[5]*B[4]+A[4]*B[5];
	C[7] = A[5]*B[3]+A[4]*B[4]+A[3]*B[5];
	C[6] = A[5]*B[2]+A[4]*B[3]+A[3]*B[4]+A[2]*B[5];
	C[5] = A[5]*B[1]+A[4]*B[2]+A[3]*B[3]+A[2]*B[4]+A[1]*B[5];
	C[4] = A[4]*B[1]+A[3]*B[2]+A[2]*B[3]+A[1]*B[4];
	C[3] = A[3]*B[1]+A[2]*B[2]+A[1]*B[3];
	C[2] = A[2]*B[1]+A[1]*B[2];
	C[1] = A[1]*B[1];
	return;
}

// polyProduct4sq {{{1
/* C becomes the polyProduct4sq of A*A */
void polyProduct4sq(const utility::vector1<Real>& A, utility::vector1<Real>& C) {
	C.resize(9);
	C[9] = A[5]*A[5];
	C[8] =             2* A[5]*A[4];
	C[7] = A[4]*A[4] + 2* A[5]*A[3];
	C[6] =             2*(A[5]*A[2]+A[4]*A[3]);
	C[5] = A[3]*A[3] + 2*(A[5]*A[1]+A[4]*A[2]);
	C[4] =             2*(A[4]*A[1]+A[3]*A[2]);
	C[3] = A[2]*A[2] + 2* A[3]*A[1];
	C[2] =             2* A[2]*A[1];
	C[1] = A[1]*A[1];
	return;
}

// polyProduct6x6 {{{1
/* C becomes the polyProduct6x6 of A and B */
void polyProduct6x6(const utility::vector1<Real>& A, const utility::vector1<Real>& B, utility::vector1<Real>& C) {
	C.resize(13);
	C[13] = A[7]*B[7];
	C[12] = A[7]*B[6]+A[6]*B[7];
	C[11] = A[7]*B[5]+A[6]*B[6]+A[5]*B[7];
	C[10] = A[7]*B[4]+A[6]*B[5]+A[5]*B[6]+A[4]*B[7];
	C[9]  = A[7]*B[3]+A[6]*B[4]+A[5]*B[5]+A[4]*B[6]+A[3]*B[7];
	C[8]  = A[7]*B[2]+A[6]*B[3]+A[5]*B[4]+A[4]*B[5]+A[3]*B[6]+A[2]*B[7];
	C[7]  = A[7]*B[1]+A[6]*B[2]+A[5]*B[3]+A[4]*B[4]+A[3]*B[5]+A[2]*B[6]+A[1]*B[7];
	C[6]  = A[6]*B[1]+A[5]*B[2]+A[4]*B[3]+A[3]*B[4]+A[2]*B[5]+A[1]*B[6];
	C[5]  = A[5]*B[1]+A[4]*B[2]+A[3]*B[3]+A[2]*B[4]+A[1]*B[5];
	C[4]  = A[4]*B[1]+A[3]*B[2]+A[2]*B[3]+A[1]*B[4];
	C[3]  = A[3]*B[1]+A[2]*B[2]+A[1]*B[3];
	C[2]  = A[2]*B[1]+A[1]*B[2];
	C[1]  = A[1]*B[1];
	return;
}

// polyProduct12x4 {{{1
/* C becomes the polyProduct12x4 of A and B */
void polyProduct12x4(const utility::vector1<Real>& A, const utility::vector1<Real>& B, utility::vector1<Real>& C) {
	C.resize(17);
	C[17] = A[13]*B[5];
	C[16] = A[13]*B[4]+A[12]*B[5];
	C[15] = A[13]*B[3]+A[12]*B[4]+A[11]*B[5];
	C[14] = A[13]*B[2]+A[12]*B[3]+A[11]*B[4]+A[10]*B[5];
	C[13] = A[13]*B[1]+A[12]*B[2]+A[11]*B[3]+A[10]*B[4]+A[9]*B[5];
	C[12] = A[12]*B[1]+A[11]*B[2]+A[10]*B[3]+A[ 9]*B[4]+A[8]*B[5];
	C[11] = A[11]*B[1]+A[10]*B[2]+A[ 9]*B[3]+A[ 8]*B[4]+A[7]*B[5];
	C[10] = A[10]*B[1]+A[ 9]*B[2]+A[ 8]*B[3]+A[ 7]*B[4]+A[6]*B[5];
	C[ 9] = A[ 9]*B[1]+A[ 8]*B[2]+A[ 7]*B[3]+A[ 6]*B[4]+A[5]*B[5];
	C[ 8] = A[ 8]*B[1]+A[ 7]*B[2]+A[ 6]*B[3]+A[ 5]*B[4]+A[4]*B[5];
	C[ 7] = A[ 7]*B[1]+A[ 6]*B[2]+A[ 5]*B[3]+A[ 4]*B[4]+A[3]*B[5];
	C[ 6] = A[ 6]*B[1]+A[ 5]*B[2]+A[ 4]*B[3]+A[ 3]*B[4]+A[2]*B[5];
	C[ 5] = A[ 5]*B[1]+A[ 4]*B[2]+A[ 3]*B[3]+A[ 2]*B[4]+A[1]*B[5];
	C[ 4] = A[ 4]*B[1]+A[ 3]*B[2]+A[ 2]*B[3]+A[ 1]*B[4];
	C[ 3] = A[ 3]*B[1]+A[ 2]*B[2]+A[ 1]*B[3];
	C[ 2] = A[ 2]*B[1]+A[ 1]*B[2];
	C[ 1] = A[ 1]*B[1];
	return;
}

// polyProduct8x8 {{{1
/* C becomes the polyProduct8x8 of A and B */
void polyProduct8x8(const utility::vector1<Real>& A, const utility::vector1<Real>& B, utility::vector1<Real>& C) {
	C.resize(17);
	C[17] = A[9]*B[9];
	C[16] = A[9]*B[8]+A[8]*B[9];
	C[15] = A[9]*B[7]+A[8]*B[8]+A[7]*B[9];
	C[14] = A[9]*B[6]+A[8]*B[7]+A[7]*B[8]+A[6]*B[9];
	C[13] = A[9]*B[5]+A[8]*B[6]+A[7]*B[7]+A[6]*B[8]+A[5]*B[9];
	C[12] = A[9]*B[4]+A[8]*B[5]+A[7]*B[6]+A[6]*B[7]+A[5]*B[8]+A[4]*B[9];
	C[11] = A[9]*B[3]+A[8]*B[4]+A[7]*B[5]+A[6]*B[6]+A[5]*B[7]+A[4]*B[8]+A[3]*B[9];
	C[10] = A[9]*B[2]+A[8]*B[3]+A[7]*B[4]+A[6]*B[5]+A[5]*B[6]+A[4]*B[7]+A[3]*B[8]+A[2]*B[9];
	C[ 9] = A[9]*B[1]+A[8]*B[2]+A[7]*B[3]+A[6]*B[4]+A[5]*B[5]+A[4]*B[6]+A[3]*B[7]+A[2]*B[8]+A[1]*B[9];
	C[ 8] = A[8]*B[1]+A[7]*B[2]+A[6]*B[3]+A[5]*B[4]+A[4]*B[5]+A[3]*B[6]+A[2]*B[7]+A[1]*B[8];
	C[ 7] = A[7]*B[1]+A[6]*B[2]+A[5]*B[3]+A[4]*B[4]+A[3]*B[5]+A[2]*B[6]+A[1]*B[7];
	C[ 6] = A[6]*B[1]+A[5]*B[2]+A[4]*B[3]+A[3]*B[4]+A[2]*B[5]+A[1]*B[6];
	C[ 5] = A[5]*B[1]+A[4]*B[2]+A[3]*B[3]+A[2]*B[4]+A[1]*B[5];
	C[ 4] = A[4]*B[1]+A[3]*B[2]+A[2]*B[3]+A[ 1]*B[4];
	C[ 3] = A[3]*B[1]+A[2]*B[2]+A[1]*B[3];
	C[ 2] = A[2]*B[1]+A[1]*B[2];
	C[ 1] = A[1]*B[1];
	return;
}

// polyProduct8sq {{{1
/* C becomes the polyProduct8sq of A*A */
void polyProduct8sq(const utility::vector1<Real>& A, utility::vector1<Real>& C) {
	C.resize(17);
	C[17] = A[9]*A[9];
	C[16] =             2* A[9]*A[8];
	C[15] = A[8]*A[8] + 2* A[9]*A[7];
	C[14] =             2*(A[9]*A[6]+A[8]*A[7]);
	C[13] = A[7]*A[7] + 2*(A[9]*A[5]+A[8]*A[6]);
	C[12] =             2*(A[9]*A[4]+A[8]*A[5]+A[7]*A[6]);
	C[11] = A[6]*A[6] + 2*(A[9]*A[3]+A[8]*A[4]+A[7]*A[5]);
	C[10] =             2*(A[9]*A[2]+A[8]*A[3]+A[7]*A[4]+A[6]*A[5]);
	C[ 9] = A[5]*A[5] + 2*(A[9]*A[1]+A[8]*A[2]+A[7]*A[3]+A[6]*A[4]);
	C[ 8] =             2*(A[8]*A[1]+A[7]*A[2]+A[6]*A[3]+A[5]*A[4]);
	C[ 7] = A[4]*A[4] + 2*(A[7]*A[1]+A[6]*A[2]+A[5]*A[3]);
	C[ 6] =             2*(A[6]*A[1]+A[5]*A[2]+A[4]*A[3]);
	C[ 5] = A[3]*A[3] + 2*(A[5]*A[1]+A[4]*A[2]);
	C[ 4] =             2*(A[4]*A[1]+A[3]*A[2]);
	C[ 3] = A[2]*A[2] + 2* A[3]*A[1];
	C[ 2] =             2* A[2]*A[1];
	C[ 1] = A[1]*A[1];
	return;
}
// }}}1

// vectorDiff {{{1
/* C becomes the difference of A - B */
/* Assumes A and B are the same size */
void vectorDiff(const utility::vector1<Real>& A, const utility::vector1<Real>& B, utility::vector1<Real>& C) {
	C.resize(A.size());
	for ( unsigned i=1; i<=A.size(); i++ ) {
		C[i]=A[i]-B[i];
	}
	return;
}


// dixonResultant {{{1
/* R becomes the Dixon resultant of the determinant of the matrix polynomial
D(t3) := R{1} + R{2} * t3 + R{3} * t3^2
*/
void dixonResultant(const utility::vector1<utility::vector1<Real> >& A,
	const utility::vector1<utility::vector1<Real> >& B,
	const utility::vector1<utility::vector1<Real> >& C,
	const utility::vector1<utility::vector1<Real> >& D,
	utility::vector1<utility::vector1<utility::vector1<Real> > >& R) {
	R.resize(DIXON_RESULTANT_SIZE);
	for ( int i=1; i<=3; i++ ) {
		R[i].resize(DIXON_SIZE);
		//first row
		R[i][1].resize(DIXON_SIZE);
		R[i][1][1]=0.0;
		R[i][1][2]=A[1][i];
		R[i][1][3]=A[2][i];
		R[i][1][4]=A[3][i];
		R[i][1][5]=0.0;
		R[i][1][6]=B[1][i];
		R[i][1][7]=B[2][i];
		R[i][1][8]=B[3][i];
		// second row
		R[i][2].resize(DIXON_SIZE);
		R[i][2][1]=A[1][i];
		R[i][2][2]=A[2][i];
		R[i][2][3]=A[3][i];
		R[i][2][4]=0.0;
		R[i][2][5]=B[1][i];
		R[i][2][6]=B[2][i];
		R[i][2][7]=B[3][i];
		R[i][2][8]=0.0;
		// third row
		R[i][3].resize(DIXON_SIZE);
		R[i][3][1]=0.0;
		R[i][3][2]=B[1][i];
		R[i][3][3]=B[2][i];
		R[i][3][4]=B[3][i];
		R[i][3][5]=0.0;
		R[i][3][6]=C[1][i];
		R[i][3][7]=C[2][i];
		R[i][3][8]=C[3][i];
		// forth row
		R[i][4].resize(DIXON_SIZE);
		R[i][4][1]=B[1][i];
		R[i][4][2]=B[2][i];
		R[i][4][3]=B[3][i];
		R[i][4][4]=0.0;
		R[i][4][5]=C[1][i];
		R[i][4][6]=C[2][i];
		R[i][4][7]=C[3][i];
		R[i][4][8]=0.0;
		// fifth row
		R[i][5].resize(DIXON_SIZE);
		R[i][5][1]=0.0;
		R[i][5][2]=0.0;
		R[i][5][3]=0.0;
		R[i][5][4]=0.0;
		R[i][5][5]=0.0;
		R[i][5][6]=D[1][i];
		R[i][5][7]=D[2][i];
		R[i][5][8]=D[3][i];
		// sixth row
		R[i][6].resize(DIXON_SIZE);
		R[i][6][1]=0.0;
		R[i][6][2]=0.0;
		R[i][6][3]=0.0;
		R[i][6][4]=0.0;
		R[i][6][5]=D[1][i];
		R[i][6][6]=D[2][i];
		R[i][6][7]=D[3][i];
		R[i][6][8]=0.0;
		// seventh row
		R[i][7].resize(DIXON_SIZE);
		R[i][7][1]=0.0;
		R[i][7][2]=D[1][i];
		R[i][7][3]=D[2][i];
		R[i][7][4]=D[3][i];
		R[i][7][5]=0.0;
		R[i][7][6]=0.0;
		R[i][7][7]=0.0;
		R[i][7][8]=0.0;
		// eighth row
		R[i][8].resize(DIXON_SIZE);
		R[i][8][1]=D[1][i];
		R[i][8][2]=D[2][i];
		R[i][8][3]=D[3][i];
		R[i][8][4]=0.0;
		R[i][8][5]=0.0;
		R[i][8][6]=0.0;
		R[i][8][7]=0.0;
		R[i][8][8]=0.0;
	}
}
// }}}1

// build_dixon_matrices {{{1
void build_dixon_matrices(
	PseudoMatrix const & A, PseudoMatrix const & B,
	PseudoMatrix const & C, PseudoMatrix const & D,
	Matrix8 & R0, Matrix8 & R1, Matrix8 & R2) {

	// The variable names used here mostly correspond to those used by:
	// Coutsias, Seok, Wester, Dill; Int. J. Quant. Chem. 106:176-189 (2006).

	Matrix8 * R[3] = {&R0, &R1, &R2};

	for ( int i = 1; i <= 3; i++ ) {
		Matrix8 & Ri = *R[i-1];
		Ri.fill(0);

		Ri(0, 1) = A[1][i];   Ri(0, 2) = A[2][i];   Ri(0, 3) = A[3][i];
		Ri(1, 0) = A[1][i];   Ri(1, 1) = A[2][i];   Ri(1, 2) = A[3][i];

		Ri(0, 5) = B[1][i];   Ri(0, 6) = B[2][i];   Ri(0, 7) = B[3][i];
		Ri(1, 4) = B[1][i];   Ri(1, 5) = B[2][i];   Ri(1, 6) = B[3][i];
		Ri(2, 1) = B[1][i];   Ri(2, 2) = B[2][i];   Ri(2, 3) = B[3][i];
		Ri(3, 0) = B[1][i];   Ri(3, 1) = B[2][i];   Ri(3, 2) = B[3][i];

		Ri(2, 5) = C[1][i];   Ri(2, 6) = C[2][i];   Ri(2, 7) = C[3][i];
		Ri(3, 4) = C[1][i];   Ri(3, 5) = C[2][i];   Ri(3, 6) = C[3][i];

		Ri(4, 5) = D[1][i];   Ri(4, 6) = D[2][i];   Ri(4, 7) = D[3][i];
		Ri(5, 4) = D[1][i];   Ri(5, 5) = D[2][i];   Ri(5, 6) = D[3][i];
		Ri(6, 1) = D[1][i];   Ri(6, 2) = D[2][i];   Ri(6, 3) = D[3][i];
		Ri(7, 0) = D[1][i];   Ri(7, 1) = D[2][i];   Ri(7, 2) = D[3][i];
	}
}

// build_sin_and_cos {{{1
void build_sin_and_cos( utility::vector1< PseudoVector > const & u, utility::vector1< PseudoVector > & sin, utility::vector1< PseudoVector >  & cos) {

	// The variable names used here mostly correspond to those used by:
	// Coutsias, Seok, Wester, Dill; Int. J. Quant. Chem. 106:176-189 (2006).

	int num_solutions = u.size();
	double u_squared, u_squared_plus_1;

	cos.resize(num_solutions);
	sin.resize(num_solutions);

	for ( int i = 1; i <= num_solutions; i++ ) {

		for ( int j = 1; j <= 3; j++ ) {
			u_squared = u[i][j] * u[i][j];
			u_squared_plus_1 = u_squared + 1;

			cos[i][j] = (1 - u_squared) / u_squared_plus_1;
			sin[i][j] = 2 * u[i][j] / u_squared_plus_1;
		}
	}
}

void dixon_eig( // {{{1
	PseudoMatrix const & A,
	PseudoMatrix const & B,
	PseudoMatrix const & C,
	PseudoMatrix const & D,
	vector1<int> const & /*order*/,
	utility::vector1< utility::fixedsizearray1<Real, 3> > & cos,
	utility::vector1< utility::fixedsizearray1<Real, 3> > & sin,
	utility::vector1< utility::fixedsizearray1<Real, 3> > & u,
	int & num_solutions) {

	// The variable names used here mostly correspond to those used by:
	// Coutsias, Seok, Wester, Dill; Int. J. Quant. Chem. 106:176-189 (2006).
	//
	// An exception is made for the script A and B variables, which represent the
	// 16x16 matrices used to setup to generalized eigenvalue problem, because
	// the non-script A and B variables are already in use.  Instead, script A
	// and B will be named U and V, respectively.

	// Fill in the Dixon matrices.
	Matrix8 R0, R1, R2;
	build_dixon_matrices(A, B, C, D, R0, R1, R2);

	// Setup the eigenvalue problem.
	Matrix16 U = Matrix16::Zero();
	Matrix16 V = Matrix16::Zero();
	Matrix8 I = Matrix8::Identity();

	U.topRightCorner(8, 8) = I;
	U.bottomLeftCorner(8, 8) = -R0;
	U.bottomRightCorner(8, 8) = -R1;

	V.topLeftCorner(8, 8) = I;
	V.bottomRightCorner(8, 8) = R2;

	// Solve the eigenvalue problem.
	SolverType solver(U, V);
	SolverType::RealEigenvalueType eigenvalues = solver.real_eigenvalues();
	SolverType::RealEigenvectorType eigenvectors = solver.real_eigenvectors();

	// Extract the closure variables.
	num_solutions = solver.num_real_solutions();
	u.resize(num_solutions);

	for ( int i = 0; i < num_solutions; i++ ) {

		if ( std::abs(eigenvectors(0, i)) > 1e-8 ) {
			u[i+1][1] = eigenvectors(1, i) / eigenvectors(0, i);
			u[i+1][2] = eigenvectors(4, i) / eigenvectors(0, i);
		} else {
			u[i+1][1] = eigenvectors(3, i) / eigenvectors(2, i);
			u[i+1][2] = eigenvectors(6, i) / eigenvectors(2, i);
		}

		u[i+1][3] = eigenvalues(i);
	}

	build_sin_and_cos(u, sin, cos);
}

void dixon_sturm( // {{{1
	const utility::fixedsizearray1<utility::fixedsizearray1<Real,3>,3 >& A,
	const utility::fixedsizearray1<utility::fixedsizearray1<Real,3>,3 >& B,
	const utility::fixedsizearray1<utility::fixedsizearray1<Real,3>,3 >& C,
	const utility::fixedsizearray1<utility::fixedsizearray1<Real,3>,3 >& D,
	const utility::vector1<int>& order,
	utility::vector1<utility::fixedsizearray1<Real,3> >& cos,
	utility::vector1<utility::fixedsizearray1<Real,3> >& sin,
	utility::vector1<utility::fixedsizearray1<Real,3> >& tau, int& nsol) {

	using utility::vector1;
	using namespace numeric::kinematic_closure;

	utility::vector1<Real> ddet (17); // Dixon determinant
	utility::vector1<Real> A0D1, A1D0, A1D2, A2D1, A0D2, A2D0, B0D1, B0D2, B1D0, B1D2, B2D0, B2D1, C0D1, C0D2, C1D0, C1D2, C2D0, C2D1, D0D2; // 1. Products of 2 quadratics
	utility::vector1<Real> Det11, Det12, Det13, Det21, Det22, Det23, Det31, Det32, Det33; // 2. Sums of quartics;
	utility::vector1<Real> D1131, D1132, D1231, D1232, D1233, D1332, D1333, D2122, D2223; // 3. Products of 2 quartics
	utility::vector1<Real> D2121, D2222, D2323; // 4. Squares of quartics
	utility::vector1<Real> Det1256, Det3478, Det1257, Det3468, Det1356, Det2478, Det1357; // 5. Sum of octics
	utility::vector1<Real> C0Det23, C1Det22, C2Det21, A0Det23, A1Det22, A2Det21; // 6. Products of a quadratic with a quartic
	utility::vector1<Real> Det678, Det123; // 7. Sums of sextics
	utility::vector1<Real> Det123_678; // 8. Product of sextics
	utility::vector1<Real> D1235_4678; // 9. Product of a quartic with a degree-12 poly
	utility::vector1<Real> D1256_3478, D1257_3468, D1356_2478; // 10. Products of octics
	utility::vector1<Real> D1357_2468; // 11. square of octic
	utility::vector1<Real> pd1, pd2, pd3; // point_value2 results
	utility::vector1<Real> pdet21, pdet22, pdet23, pdet31, pdet32, pdet33, pd0d2; // point_value4 results
	utility::vector1<Real> cram12, crn112, crn262; // point_value6 results
	utility::vector1<Real> cram22, cram32, cram42, cram52, cram122, cram132, cram222, cram232, cram242, cram252; // point_value8 results
	utility::vector1<Real> cram11, cram21, cram31, cram41, cram51, crden; // sums and differences of products of point_value results used to find denominator crden
	utility::vector1<Real> crn111, crn121, crn131, crn122, crn132; // sums and differences of products of point_value results used to find Cramer numerator of tau[1]
	utility::vector1<Real> crn221, crn222, crn231, crn232, crn241, crn242, crn251, crn252, crn261; // products of point_value results used to find Cramer numerator of tau[2]
	utility::vector1<Real> cram_num1, cram_num2; // Cramer numerators for tau[1] and tau[2]
	utility::vector1<Real> z1, z2; // tau[1] and tau[2] numerator divided element-wise by denominator crden
	//utility::vector1<utility::vector1<Real> > tau (16); // half tangents
	utility::vector1<utility::vector1<Real> > tsq (16); // tau.^2
	utility::vector1<utility::vector1<Real> > tsq1 (16); // tsq+1
	// sturm initialization constants
	Real tol_secant=1.0e-15;
	int max_iter_sturm=100, max_iter_secant=20;
	int p_order=16; // order of polynomial for sturm
	utility::vector1<Real> Q (17); // polynomial coefficients for sturm solver
	utility::vector1<Real> roots; // roots found by sturm solver

	////----------------------------------------------------------------------------
	// A: Characteristic polynomial via Lagrange 4x4-expansion of Dixon determinant.
	//    Roots of characteristic polynomial will give tau[3]. In B, we'll use
	//    Cramer's rule to derive tau[1] and tau[2] from tau[3].
	////----------------------------------------------------------------------------

	//-----------------------------------
	//A1. Products of 2 quadratics
	//   (19 * 13ops = 247ops]
	//-----------------------------------
	polyProduct2x2(A[1],D[2],A0D1); // A0*D1
	polyProduct2x2(A[2],D[1],A1D0); // A1*D0:
	polyProduct2x2(A[2],D[3],A1D2); // A1*D2:
	polyProduct2x2(A[3],D[2],A2D1); // A2*D1:
	polyProduct2x2(A[1],D[3],A0D2); // A0*D2:
	polyProduct2x2(A[3],D[1],A2D0); // A2*D0:
	polyProduct2x2(B[1],D[2],B0D1); // B0*D1:
	polyProduct2x2(B[1],D[3],B0D2); // B0*D2:
	polyProduct2x2(B[2],D[1],B1D0); // B1*D0:
	polyProduct2x2(B[2],D[3],B1D2); // B1*D2:
	polyProduct2x2(B[3],D[1],B2D0); // B2*D0:
	polyProduct2x2(B[3],D[2],B2D1); // B2*D1:
	polyProduct2x2(C[1],D[2],C0D1); // C0*D1:
	polyProduct2x2(C[1],D[3],C0D2); // C0*D2:
	polyProduct2x2(C[2],D[1],C1D0); // C1*D0:
	polyProduct2x2(C[2],D[3],C1D2); // C1*D2:
	polyProduct2x2(C[3],D[1],C2D0); // C2*D0:
	polyProduct2x2(C[3],D[2],C2D1); // C2*D1:
	polyProduct2x2(D[1],D[3],D0D2); // D0*D2:

	//-----------------------------------
	//A2. Sums of  quartics
	//   (9 * 5ops = 45ops)
	//-----------------------------------

	vectorDiff(A0D1, A1D0, Det11);
	vectorDiff(A0D2, A2D0, Det12);
	vectorDiff(A1D2, A2D1, Det13);
	vectorDiff(B0D1, B1D0, Det21);
	vectorDiff(B0D2, B2D0, Det22);
	vectorDiff(B1D2, B2D1, Det23);
	vectorDiff(C0D1, C1D0, Det31);
	vectorDiff(C0D2, C2D0, Det32);
	vectorDiff(C1D2, C2D1, Det33);

	//-----------------------------------
	//A3. Products of 2 quartics
	//   (9 * 41ops = 369ops)
	//-----------------------------------

	polyProduct4x4(Det11,Det31,D1131); // Det11*Det31:
	polyProduct4x4(Det11,Det32,D1132); // Det11*Det32:
	polyProduct4x4(Det12,Det31,D1231); // Det12*Det31:
	polyProduct4x4(Det12,Det32,D1232); // Det12*Det32:
	polyProduct4x4(Det12,Det33,D1233); // Det12*Det33:
	polyProduct4x4(Det13,Det32,D1332); // Det13*Det32:
	polyProduct4x4(Det13,Det33,D1333); // Det13*Det33:
	polyProduct4x4(Det21,Det22,D2122); // Det21*Det22:
	polyProduct4x4(Det22,Det23,D2223); // Det22*Det23:

	//-----------------------------------
	//A4. Squares of quartics
	//   (3 * 28ops = 84ops)
	//-----------------------------------

	polyProduct4sq(Det21, D2121); // Det21^2:
	polyProduct4sq(Det22, D2222); // Det22^2:
	polyProduct4sq(Det23, D2323); // Det23^2:

	//-----------------------------------
	//A5. sums of octics
	//   (7 * 9ops = 63ops)
	//-----------------------------------
	// second factor of 5th product identical with first
	// not used (shown for clarity)
	//   Det2468 =  Det1357;

	vectorDiff(D2121,D1131,Det1256);
	vectorDiff(D2323,D1333,Det3478);
	vectorDiff(D2122,D1132,Det1257);
	vectorDiff(D2223,D1332,Det3468);
	vectorDiff(D2122,D1231,Det1356);
	vectorDiff(D2223,D1233,Det2478);
	vectorDiff(D2222,D1232,Det1357);

	//-----------------------------------
	//A6. products of a quadratic with a quartic
	//   (6 * 23 = 138ops)
	//-----------------------------------

	polyProduct4x2(Det23,C[1],C0Det23); // C0*Det23:
	polyProduct4x2(Det22,C[2],C1Det22); // C1*Det22:
	polyProduct4x2(Det21,C[3],C2Det21); // C2*Det21:
	polyProduct4x2(Det23,A[1],A0Det23); // A0*Det23:
	polyProduct4x2(Det22,A[2],A1Det22); // A1*Det22:
	polyProduct4x2(Det21,A[3],A2Det21); // A2*Det21:

	//-----------------------------------
	//A7. sums of sextics
	//   (4 * 7ops = 28ops)
	//-----------------------------------

	Det678.resize(PP4x2_VECSIZE); // same size as polyProduct4x2 output
	Det123.resize(PP4x2_VECSIZE); // same size as polyProduct4x2 output
	// compute Det678 and Det123
	for ( int i=1; i<=PP4x2_VECSIZE; i++ ) {
		Det678[i]=(-C0Det23[i] + C1Det22[i] - C2Det21[i]);
		Det123[i]=( A0Det23[i] - A1Det22[i] + A2Det21[i]);
	}

	//-----------------------------------
	//A8. product of sextics
	//   ( 85ops)
	//-----------------------------------

	polyProduct6x6(Det123,Det678,Det123_678); // Det123*Det678;

	//-----------------------------------
	//A9. product of a quartic with a degree-12 poly
	//   ( 113ops)
	//-----------------------------------

	polyProduct12x4(Det123_678,D0D2,D1235_4678); // D0D2*Det123_678:

	//-----------------------------------
	//A10. products of octics
	//   ( 3 * 145ops = 435ops)
	//-----------------------------------

	polyProduct8x8(Det1256,Det3478,D1256_3478); // Det1256*Det3478:
	polyProduct8x8(Det1257,Det3468,D1257_3468); // Det1257*Det3468:
	polyProduct8x8(Det1356,Det2478,D1356_2478); // Det1356*Det2478:

	//-----------------------------------
	//A11. square of octic -- since Det2468 = Det1357
	//    ( 88ops)
	//-----------------------------------
	// 6th product is identical to the first
	// Not used (shown for clarity)
	// D1567_2348 = D1235_4678;

	polyProduct8sq(Det1357,D1357_2468); // Det1357*Det2468:

	//-----------------------------------
	//A12.  Sums of 16th degree polynomials
	//     ( 5 * 17ops = 85ops)
	//-----------------------------------

	for ( unsigned i=1; i<=ddet.size(); i++ ) {
		ddet[i] = (-2*D1235_4678[i] + D1256_3478[i] - D1257_3468[i] - D1356_2478[i] + D1357_2468[i]);
	}

	//-----------------------------------
	//A13.  Solve for the roots
	//-----------------------------------

	// Get the coefficients
	for ( unsigned i=1; i<=ddet.size(); i++ ) {
		Q[i]=ddet[i] / ddet[17];
	}

	// Get the roots and put them into tau[3]
	initialize_sturm(&tol_secant, &max_iter_sturm, &max_iter_secant);
	solve_sturm(p_order, nsol, Q, roots);
	if ( nsol==0 ) return; // no solns found

	// remove negative roots // DJM: 7-10-2007: KEEP ALL ROOTS!
	//utility::vector1<Real>::iterator rootIterator=roots.end();
	//while(rootIterator != roots.begin()) {
	// rootIterator--;
	// if (*rootIterator < 0) {
	//  roots.erase(rootIterator);
	// }
	//}

	nsol=roots.size();
	tau.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		tau[i][order[3]]=roots[i];  // restoring column-major indexing!
	}

	////---------------------------------------------------------------------
	// B: Determination of tau[1] and tau[2] from tau[3] using Cramers' rule
	////---------------------------------------------------------------------

	//-----------------------------------
	//B1.  Compute the denominator
	//-----------------------------------

	point_value4(Det23,roots,pdet23);
	point_value4(D0D2,roots,pd0d2);
	cram11.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram11[i] = pdet23[i]*pd0d2[i];
	}
	point_value6(Det678, roots, cram12);
	point_value4(Det31, roots, pdet31);
	point_value2(D[2], roots, pd2);
	cram21.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram21[i] = -pdet31[i]*pd2[i];
	}
	point_value8(Det3478, roots, cram22);
	point_value4(Det32, roots, pdet32);
	cram31.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram31[i] = -pdet32[i]*pd2[i];
	}
	point_value8(Det3468, roots, cram32);
	point_value2(D[3], roots, pd3);
	cram41.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram41[i] = -pdet31[i]*pd3[i];
	}
	point_value8(Det2478, roots, cram42);
	cram51.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram51[i] = -pdet32[i]*pd3[i];
	}
	point_value8(Det1357, roots, cram52);
	crden.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crden[i] = -cram11[i]*cram12[i] + cram21[i]*cram22[i] - cram31[i]*cram32[i] - cram41[i]*cram42[i] + cram51[i]*cram52[i];
	}

	//-------------------------------------------------------------
	//B2.  Compute the numerator for tau[1] and store tau[1] values
	//-------------------------------------------------------------

	point_value4(Det22, roots, pdet22);
	crn111.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn111[i] = pdet22[i]*pd0d2[i];
	}
	point_value6(Det678, roots, crn112);
	point_value2(D[1], roots, pd1);
	crn121.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn121[i] = -pdet31[i]*pd1[i];
	}
	point_value8(Det3478, roots, crn122);
	crn131.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn131[i] = -pdet32[i]*pd1[i];
	}
	point_value8(Det3468, roots, crn132);
	cram_num1.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram_num1[i] = -crn111[i]*crn112[i]+crn121[i]*crn122[i]-crn131[i]*crn132[i];
		tau[i][order[1]] = -cram_num1[i] / crden[i]; // restoring column-major indexing!
	}

	//-------------------------------------------------------------
	//B3.  Compute the numerator for tau[2] and store tau[2] values
	//-------------------------------------------------------------

	point_value4(Det21, roots, pdet21);
	crn221.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn221[i] = -pdet21[i]*pd2[i];
	}
	point_value8(Det3478, roots, crn222);
	crn231.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn231[i] = -pdet21[i]*pd3[i];
	}
	point_value8(Det3468, roots, crn232);
	crn241.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn241[i] = -pdet22[i]*pd2[i];
	}
	point_value8(Det2478, roots, crn242);
	crn251.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn251[i] = -pdet22[i]*pd3[i];
	}
	point_value8(Det1357, roots, crn252);
	point_value4(Det33, roots, pdet33);
	crn261.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		crn261[i] = -pdet33[i]*pd0d2[i];
	}
	point_value6(Det123, roots, crn262);
	cram_num2.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		cram_num2[i] = crn221[i]*crn222[i] - crn231[i]*crn232[i] - crn241[i]*crn242[i] + crn251[i]*crn252[i] - crn261[i]*crn262[i];
		tau[i][order[2]] = -cram_num2[i] / crden[i]; // DJM: division inside for-loop!
	}                                              // and restoring column-major indexing!

	//----------------------------------------
	//B4. Compute the output sines and cosines
	//----------------------------------------
	tsq.resize(nsol);
	tsq1.resize(nsol);
	cos.resize(nsol);
	sin.resize(nsol);
	for ( int i=1; i<=nsol; i++ ) {
		tsq[i].resize(3);
		tsq1[i].resize(3);
		for ( int j=1; j<=3; j++ ) {
			tsq[i][j]=(tau[i][j])*(tau[i][j]);
			tsq1[i][j]=tsq[i][j]+1;
			cos[i][j]=(1-tsq[i][j]) / tsq1[i][j];  // DJM: division inside
			sin[i][j]=2 * tau[i][j] / tsq1[i][j];  // of for-loops!
		}
	}

	return;
}
// }}}1

// test_point_value2 {{{1
void test_point_value2() {
	utility::vector1<Real> t(16), C;
	utility::fixedsizearray1< Real, 3 > A;
	A[1]=0.0;
	A[2]=-60;
	A[3]=0.0;
	t[1]=7.4248;
	t[2]=5.3100;
	t[3]=0.0;
	t[4]=0.0;
	t[5]=0.0;
	t[6]=0.0;
	t[7]=0.0;
	t[8]=0.0;
	t[9]=0.0;
	t[10]=0.0;
	t[11]=0.0;
	t[12]=0.0;
	t[13]=0.0;
	t[14]=0.0;
	t[15]=0.0;
	t[16]=0.0;
	point_value2(A,t,C);
	printVector(C);
	// output should be { -445.4891, -318.6002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0, 0.0, 0.0}
}

// test_polyProduct6x6 {{{1
void test_polyProduct6x6() {
	utility::vector1<Real> A (7), B (7), C;
	A[1]=366130000.0;
	A[2]=0.0;
	A[3]=316600000.0;
	A[4]=0.0;
	A[5]=-47950000.0;
	A[6]=0.0;
	A[7]=920000.0;
	B[1]=1790200000.0;
	B[2]=0.0;
	B[3]=2113200000.0;
	B[4]=0.0;
	B[5]=317100000.0;
	B[6]=0.0;
	B[7]=-8300000.0;
	polyProduct6x6(A,B,C);
	printVector(C);
	// output should be 1.0e+18*{0.6554, 0.0, 1.3405, 0.0, 0.6993, 0.0, -0.0023, 0.0, -0.0159, 0.0}, 6.89717e+14, 0.0, -7.636e+12
	return;
}

// test_polyProduct4sq {{{1
void test_polyProduct4sq() {
	utility::vector1<Real> A (5), C;
	A[1]=126600.0;
	A[2]=0.0;
	A[3]=-28200.0;
	A[4]=0.0;
	A[5]=-2890.0;
	polyProduct4sq(A,C);
	printVector(C);
	// output should be 1.0e+10*{1.6027, 0.0, -0.7140, 0.0, 0.0063, 0.0, 0.0163, 0.0, 0.0008}
	return;
}

// test_polyProduct4x4 {{{1
void test_polyProduct4x4() {
	utility::vector1<Real> A (5), B (5), C;
	A[1]=126940.0;
	A[2]=0.0;
	A[3]=-28110.0;
	A[4]=0.0;
	A[5]=-2890.0;
	B[1]=0.0;
	B[2]=172730.0;
	B[3]=0.0;
	B[4]=20800.0;
	B[5]=0.0;
	polyProduct4x4(A,B,C);
	printVector(C);
	// output should be 1.0e+10*{0.0, 2.1927, 0.0, -0.2216, 0.0, -0.1085, 0.0, -0.0060, 0.0}
	return;
}

// test_polyProduct4x2 {{{1
void test_polyProduct4x2() {
	utility::vector1<Real> A (5), C;
	utility::fixedsizearray1< Real, 3 > B;
	A[1]=0.0;
	A[2]=-108950.0;
	A[3]=0.0;
	A[4]=43180.0;
	A[5]=0.0;
	B[1]=0.0;
	B[2]=-1799.10;
	B[3]=0.0;
	polyProduct4x2(A,B,C);
	printVector(C);
	// output should be 1.0e+08 * {0.0, 0.0, 1.9601, 0.0, -0.7769, 0.0, 0.0}
	return;
}

// test_polyProduct2x2 {{{1
void test_polyProduct2x2() {
	utility::fixedsizearray1<Real,3> A, B;
	utility::vector1< Real > C;
	A[1]=15.999410431075338;
	A[2]=0.0;
	A[3]=-8.999904923464781;
	B[1]=72.015965384650187;
	B[2]=0.0;
	B[3]=7.000641402062060;
	polyProduct2x2(A,B,C);
	printVector(C);
	// output should be 1.03e+03 * {1.152212987779133, 0.0, -0.536130706361013, 0.0, -0.063005107021830}
	return;
}

// test_dixon() {{{1
void test_dixon() {
	Real Avals[]={-0.061472096577788, -2.424098126176366, -0.016332287075674, 2.421954432062721, -0.090242315482933, 0.643479843225995, 0.061472096577788, -0.646358402699405, 0.016332287075674};
	Real Bvals[]={-0.000625561301665, 0.0, 1.070810990722253,  0.054344479784912, 0.0, 0.054308028868740, -1.071191403371629, 0.0, 0.000963218011333};
	Real Cvals[]={0.063394336809916, 2.498922707477288, -0.025585174231195, -2.497689253369826, 0.177963037883907, 1.008036647097698, -0.063394336809916, -1.006882302775246, 0.025585174231195};
	Real Dvals[]={0.673457026339820, -0.015821144672541, 1.005285612756159, -0.078123031130501, 3.277064399971811, 0.078123031130501, 1.004894243686905, 0.075272539656448, -0.573852487168651};
	//utility::vector1<utility::vector1<utility::vector1<Real> > > P(3); // dixon input
	utility::vector1<utility::vector1<Real> > A (3);
	utility::vector1<utility::vector1<Real> > B (3);
	utility::vector1<utility::vector1<Real> > C (3);
	utility::vector1<utility::vector1<Real> > D (3);
	utility::vector1<int> order (3);  // dixon input
	utility::vector1<utility::vector1<Real> > tau; // dixon output
	utility::vector1<utility::vector1<Real> > cos; // dixon output
	utility::vector1<utility::vector1<Real> > sin; // dixon output
	int nsol = 0;

	// Allocate A,B,C,D
	for ( int i=1; i<=3; i++ ) {
		A[i].resize(3);
		B[i].resize(3);
		C[i].resize(3);
		D[i].resize(3);
	}

	// Fill A,B,C,D
	int n=0;
	for ( int i=1; i<=3; i++ ) {
		for ( int j=1; j<=3; j++ ) {
			A[i][j]=Avals[n];
			B[i][j]=Bvals[n];
			C[i][j]=Cvals[n];
			D[i][j]=Dvals[n];
			n++;
		}
	}

	// initialize order
	order[1]=1;
	order[2]=2;
	order[3]=3;

	// compute Dixon resultant
	//dixon(A, B, C, D, order, cos, sin, tau, nsol);

	cout << "number of solutions" << nsol << std::endl;
	printMatrix(cos);
	printMatrix(sin);
	return;
}
// }}}1

// main (commented out) {{{1
// int main(int argc, char** argv)
// {
//   //test_polyProduct2x2();
//   //test_polyProduct4x2();
//   //test_polyProduct4x4();
//   //test_polyProduct4sq();
//   //test_polyProduct6x6();
//  //test_point_value2();
//  // for (int i=1; i<10000; i++) {
//  test_dixon();
//  //}
// }
// }}}1

} // end namespace kinematic_closure
} // end namespace numeric
