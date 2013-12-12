// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   DixonResultantTests.cxxtest.hh
/// @brief  Unit tests for the Dixon resultant.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Headers {{{1

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <numeric/types.hh>
#include <numeric/kinematic_closure/dixon.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/linear_algebra/GeneralizedEigenSolver.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

// C++ headers
#include <iostream>
#include <iomanip>
#include <string>
#include <unistd.h>

// Defines {{{1
#define setup setUp

// Namespaces {{{1
using namespace std;
using namespace core;
using namespace numeric::kinematic_closure;

using utility::vector1;
using numeric::Real;
using numeric::linear_algebra::GeneralizedEigenSolver;
using Eigen::Matrix;
using Eigen::Dynamic;

// Typedefs {{{1
typedef Matrix<double, 8, 8> Matrix8;
typedef Matrix<double, 16, 16> Matrix16;
typedef Matrix<double, Dynamic, 1, 0, 16, 1> Vector16;
typedef GeneralizedEigenSolver<Matrix16> SolverType;

typedef vector1<Real> PseudoVector;
typedef vector1<PseudoVector> PseudoMatrix;
// }}}1

class DixonResultantTests : public CxxTest::TestSuite {
public:

	void setup() { // {{{1
		core_init();
    chdir("../../test/numeric/kinematic_closure");
	}
	// }}}1

	void dont_test_given_R0_R1_R2() { // {{{1
		Matrix8 R0, R1, R2;
		Matrix8 I = Matrix8::Identity();

		Eigen::IOFormat scipy(8, 0, ", ", ",\n", "[", "]", "[", "]");
		string skipl = "\n\n";

		R0 <<
									 0,     1.74673017409,                 0,   -0.385188709827,                 0,                 0,                 0,    0.886266959645,
			 1.74673017409,    0.646740198748,   -0.385188709827,    -2.74541524944,                 0,                 0,    0.886266959645,    0.102021286235,
			0.646740198748,    -1.53905498197,    -2.74541524944,   -0.858501270662,                 0,                 0,    0.102021286235,    0.746600095747,
			-1.53905498197,                 0,   -0.858501270662,                 0,                 0,                 0,    0.746600095747,                 0,
									 0,   -0.385188709827,                 0,    -1.51390472194,                 0,    0.886266959645,                 0,                 0,
		 -0.385188709827,    -2.74541524944,    -1.51390472194,   -0.343410160552,    0.886266959645,    0.102021286235,                 0,                 0,
			-2.74541524944,   -0.858501270662,   -0.343410160552,     1.43383899338,    0.102021286235,    0.746600095747,                 0,                 0,
		 -0.858501270662,                 0,     1.43383899338,                 0,    0.746600095747,                 0,                 0,                 0;

		R1 <<
									 0,    0.301512618879,                 0,   -0.350315204978,                 0,                 0,                 0,  -0.0543075131269,
			0.301512618879,    0.388220669491,   -0.350315204978,    0.190485294838,                 0,                 0,  -0.0543075131269,     3.53706323054,
			0.388220669491,     1.02615980037,    0.190485294838,   -0.317475333953,                 0,                 0,     3.53706323054,  -0.0515545052508,
			 1.02615980037,                 0,   -0.317475333953,                 0,                 0,                 0,  -0.0515545052508,                 0,
									 0,   -0.350315204978,                 0,   -0.663485945612,                 0,  -0.0543075131269,                 0,                 0,
		 -0.350315204978,    0.190485294838,   -0.663485945612,    0.912942501326,  -0.0543075131269,     3.53706323054,                 0,                 0,
			0.190485294838,   -0.317475333953,    0.912942501326,    0.151623847683,     3.53706323054,  -0.0515545052508,                 0,                 0,
		 -0.317475333953,                 0,    0.151623847683,                 0,  -0.0515545052508,                 0,                 0,                 0;


		R2 <<
									 0,   -0.800486445641,                 0,   0.0241631939177,                 0,                 0,                 0,  -0.0551498510613,
		 -0.800486445641,    -0.33761472333,   0.0241631939177,     1.00287787509,                 0,                 0,  -0.0551498510613,   -0.102021286235,
			-0.33761472333,    0.512749686414,     1.00287787509,    0.197060425169,                 0,                 0,   -0.102021286235,   -0.147093522617,
			0.512749686414,                 0,    0.197060425169,                 0,                 0,                 0,   -0.147093522617,                 0,
									 0,   0.0241631939177,                 0,    0.477903664399,                 0,  -0.0551498510613,                 0,                 0,
		 0.0241631939177,     1.00287787509,    0.477903664399,  -0.0826531058538,  -0.0551498510613,   -0.102021286235,                 0,                 0,
			 1.00287787509,    0.197060425169,  -0.0826531058538,   -0.366971525496,   -0.102021286235,   -0.147093522617,                 0,                 0,
			0.197060425169,                 0,   -0.366971525496,                 0,   -0.147093522617,                 0,                 0,                 0;

		// Note that if I forget to transpose these matrices, I still get the right 
		// eigenvalues, but I get the wrong eigenvectors.  It's very surprising to 
		// me that I would even get the eigenvalues right.  
		R0.transposeInPlace();
		R1.transposeInPlace();
		R2.transposeInPlace();

		Matrix16 A, B;

		A.fill(0);
		B.fill(0);

		//B.topLeftCorner(8, 8) = I;
		//B.bottomRightCorner(8, 8) = R2;

		//A.topRightCorner(8, 8) = I;
		//A.bottomLeftCorner(8, 8) = -R0;
		//A.bottomRightCorner(8, 8) = -R1;

		B.block<8,8>(0, 0) = I;
		B.block<8,8>(8, 8) = R2;

		A.block<8,8>(0, 8) = I;
		A.block<8,8>(8, 0) = -R0;
		A.block<8,8>(8, 8) = -R1;

		//cout << "R0 = \\" << endl << R0.format(scipy) << skipl;
		//cout << "R1 = \\" << endl << R1.format(scipy) << skipl;
		//cout << "R2 = \\" << endl << R2.format(scipy) << skipl;
		//cout << endl << endl;
		//cout << "A = \\" << endl << A.format(scipy) << skipl;
		//cout << "B = \\" << endl << B.format(scipy) << skipl;

		SolverType solver(A, B);
		SolverType::RealEigenvalueType eigenvalues;
		SolverType::RealEigenvectorType eigenvectors;

		eigenvalues = solver.real_eigenvalues();
		eigenvectors = solver.real_eigenvectors();
		int num_solutions = solver.num_real_solutions();

		Vector16 u1(num_solutions), u2(num_solutions), u3(num_solutions);

		u1.fill(0);
		u2.fill(0);
		u3.fill(0);

		cout << "------------------------------------------------------------";
		cout << endl << endl;

		cout << "eigenvectors = \\" << endl;
		cout << eigenvectors.format(scipy) << endl << endl;

		for (int i = 0; i < num_solutions; i++) {
			if (abs(eigenvectors(1, i)) > 1e-8) {
				cout << "Using: 1, 2, 5" << endl;
				u1(i) = eigenvectors(2, i) / eigenvectors(1, i);
				u2(i) = eigenvectors(5, i) / eigenvectors(1, i);
			}
			else {
				cout << "Using: 3, 4, 7" << endl;
				u1(i) = eigenvectors(4, i) / eigenvectors(3, i);
				u2(i) = eigenvectors(7, i) / eigenvectors(3, i);
			}
			u3(i) = eigenvalues(i);
		}
		cout << endl;

		cout << "u3 (via eigenvalues) = \\" << endl;
		cout << u3.format(scipy) << endl << endl;

		cout << "u2 (via eigenvectors) = \\" << endl;
		cout << u2.format(scipy) << endl << endl;

		cout << "u1 (via eigenvectors) = \\" << endl;
		cout << u1.format(scipy) << endl << endl;

	}
	// }}}1
	void dont_test_given_A_B_C_D() { // {{{1
		PseudoMatrix A, B, C, D;
		PseudoMatrix u, sin, cos;
		vector1<int> order;		// Not really used.
		int num_solutions;

		A.resize(3); B.resize(3); C.resize(3); D.resize(3);

		for (int i = 1; i <= 3; i++) {
			A[i].resize(3, 0);
			B[i].resize(3, 0);
			C[i].resize(3, 0);
			D[i].resize(3, 0);
		}

		// These inputs give an output that is nearly 180 deg.  This can provoke 
		// some numerical instabilities.

		A[1][1] =  1.74673020;   A[1][2] =  0.30151262;   A[1][3] = -0.80048645;
		A[2][1] =  0.64674020;   A[2][2] =  0.38822067;   A[2][3] = -0.33761472;
		A[3][1] = -1.53905500;   A[3][2] =  1.02615980;   A[3][3] =  0.51274969;

		B[1][1] = -0.38518871;   B[1][2] = -0.35031520;   B[1][3] =  0.02416319;
		B[2][1] = -2.74541520;   B[2][2] =  0.19048529;   B[2][3] =  1.00287790;
		B[3][1] = -0.85850127;   B[3][2] = -0.31747533;   B[3][3] =  0.19706043;

		C[1][1] = -1.51390470;   C[1][2] = -0.66348595;   C[1][3] =  0.47790366;
		C[2][1] = -0.34341016;   C[2][2] =  0.91294250;   C[2][3] = -0.08265311;
		C[3][1] =  1.43383900;   C[3][2] =  0.15162385;   C[3][3] = -0.36697153;

		D[1][1] =  0.88626696;   D[1][2] = -0.05430751;   D[1][3] = -0.05514985;
		D[2][1] =  0.10202129;   D[2][2] =  3.53706320;   D[2][3] = -0.10202129;
		D[3][1] =  0.74660010;   D[3][2] = -0.05155451;   D[3][3] = -0.14709352;

		order.resize(3);
		order[1] = 1; order[2] = 2; order[3] = 3;

		dixon_eig(A, B, C, D, order, cos, sin, u, num_solutions);

		cout << "u (accurate) = \\" << endl;
		printMatrix(u);
		cout << endl;

		dixon_sturm(A, B, C, D, order, cos, sin, u, num_solutions);

		cout << "u (inaccurate) = \\" << endl;
		printMatrix(u);
		cout << endl;

	}

	void dont_test_for_null_input() { // {{{1
		PseudoMatrix A, B, C, D;
		PseudoMatrix u, sin, cos;
		vector1<int> order;		// Not really used.
		int num_solutions;

		Eigen::IOFormat scipy(8, 0, ", ", ",\n", "[", "]", "[", "]");
		string skipl = "\n\n";

		A.resize(3); B.resize(3); C.resize(3); D.resize(3);

		for (int i = 1; i <= 3; i++) {
			A[i].resize(3, 0);
			B[i].resize(3, 0);
			C[i].resize(3, 0);
			D[i].resize(3, 0);
		}

		A[1][1] =  0.00000000;   A[1][2] =  0.00000000;   A[1][3] =  0.00000000;
		A[2][1] =  0.00000000;   A[2][2] =  0.00000000;   A[2][3] =  0.00000000;
		A[3][1] =  0.00000000;   A[3][2] =  0.00000000;   A[3][3] =  0.00000000;

		B[1][1] =  0.00000000;   B[1][2] =  0.00000000;   B[1][3] =  0.00000000;
		B[2][1] =  0.00000000;   B[2][2] =  0.00000000;   B[2][3] =  0.00000000;
		B[3][1] =  0.00000000;   B[3][2] =  0.00000000;   B[3][3] =  0.00000000;

		C[1][1] =  0.00000000;   C[1][2] =  0.00000000;   C[1][3] =  0.00000000;
		C[2][1] =  0.00000000;   C[2][2] =  0.00000000;   C[2][3] =  0.00000000;
		C[3][1] =  0.00000000;   C[3][2] =  0.00000000;   C[3][3] =  0.00000000;

		D[1][1] =  0.00000000;   D[1][2] =  0.00000000;   D[1][3] =  0.00000000;
		D[2][1] =  0.00000000;   D[2][2] =  0.00000000;   D[2][3] =  0.00000000;
		D[3][1] =  0.00000000;   D[3][2] =  0.00000000;   D[3][3] =  0.00000000;

		cout << "A = \\" << endl; printMatrix(A); cout << endl;
		cout << "B = \\" << endl; printMatrix(B); cout << endl;
		cout << "C = \\" << endl; printMatrix(C); cout << endl;
		cout << "D = \\" << endl; printMatrix(D); cout << endl;

		dixon_eig(A, B, C, D, order, cos, sin, u, num_solutions);

		cout << "u = \\" << endl; printMatrix(u); cout << endl;
	
	}
	// }}}1

};

