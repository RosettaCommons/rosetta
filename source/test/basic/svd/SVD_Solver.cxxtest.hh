// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    basic/svd/SVD_Solver.cxxtest.hh
///
/// @brief   unit test for SVD solver class
///
/// @details
/// Solve over-determined set of linear equation to minimize ||A x - b||^2, using Singular Value Decomposition (SVD) method.
/// Last modified 05/06/2016
///
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
//#include <test/core/init_util.hh>

// Unit headers
#include <basic/svd/SVD_Solver.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Objexx headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <string>

static basic::Tracer TR("basic.svd.SVD_Solver.cxxtest");

using namespace core;
using namespace basic::svd;

class SVDSolverTest : public CxxTest::TestSuite {
public:
	void setUp() {}

	void tearDown() {}

	void test_Solv_SLE_with_unique_solution() {

		// MatA and vecB as utility::vector1
		utility::vector1< utility::vector1 < Real > > MatA(4, utility::vector1<Real>(3));
		MatA[1][1] = 4;  MatA[1][2] = -2; MatA[1][3] = 1;
		MatA[2][1] = -1; MatA[2][2] = 3;  MatA[2][3] = 4;
		MatA[3][1] = 5;  MatA[3][2] = -1; MatA[3][3] = 3;
		MatA[4][1] = 2;  MatA[4][2] = 1;  MatA[4][3] = 6;

		utility::vector1< Real > vecB(4);
		vecB[1] = 15; vecB[2] = 15; vecB[3] = 26; vecB[4] = 33;

		utility::vector1< Real > vecX(3);
		Size vecB_size(4);
		Size vecX_size(3);

		// Arguments: data size, parameters size
		SVD_Solver solver(vecB_size,vecX_size);
		solver.set_matrix_A(MatA);
		solver.set_vector_b(vecB);

		// Matrix can be decomposed and SLE solved only
		// after matrix A and vector b are set.
		solver.run_decomp_svd();
		solver.run_solve_svd();
		vecX = solver.get_svd_solution();

		TS_ASSERT_DELTA(vecX[1], 2.0, 0.00001);
		TS_ASSERT_DELTA(vecX[2], -1.0, 0.00001);
		TS_ASSERT_DELTA(vecX[3], 5.0, 0.00001);

		// Return the score SQRT(||Ax-b||^2). Can be called only when SLE is solved.
		// You need the give the original matrix A as argument.
		Real score = solver.run_score_svd_on_matrix(MatA);
		TS_ASSERT_DELTA(score, 0.0, 0.00001);
	}

	void test_Solv_SLE_no_unique_solution() {
		// MatA as FArray2D and vecB as FArray1D
		ObjexxFCL::FArray2D_double MatA(4,3);
		MatA(1,1) = 4;  MatA(1,2) = -2; MatA(1,3) = 2;
		MatA(2,1) = -1; MatA(2,2) = 3;  MatA(2,3) = 4;
		MatA(3,1) = 5;  MatA(3,2) = -1; MatA(3,3) = -3;
		MatA(4,1) = 2;  MatA(4,2) = -1; MatA(4,3) = 6;

		ObjexxFCL::FArray1D_double vecB(4);
		vecB(1) = 15; vecB(2) = 15; vecB(3) = 26; vecB(4) = 33;
		utility::vector1< Real > vecX(3);
		Size vecB_size(4);
		Size vecX_size(3);

		SVD_Solver solver(vecB_size,vecX_size);
		solver.set_matrix_A(MatA);
		solver.set_vector_b(vecB);

		solver.run_decomp_svd();
		solver.run_solve_svd();
		vecX = solver.get_svd_solution();

		TS_ASSERT_DELTA(vecX[1], 6.91468, 0.00001);
		TS_ASSERT_DELTA(vecX[2], 4.43658, 0.00001);
		TS_ASSERT_DELTA(vecX[3], 2.78311, 0.00001);

		// Return the score SQRT(||Ax-b||^2). Can be called only when SLE is solved.
		// You need the give the original matrix A as argument
		Real score = solver.run_score_svd_on_matrix(MatA);
		TS_ASSERT_DELTA(score, 12.6221, 0.00001);
	}
};
