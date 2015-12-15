// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////
///
/// @file basic/svd/SVD_Solver.hh
///
/// @brief SVD solver class
///
/// @details Solve over-determined set of linear equation to minimize ||A x - b||^2, using Singular Value Decomposition (SVD) method.
///
/// @param
/// Specify the size of the problem in the constructor (M is the number of equations, N is the number of parameters to fit)
/// M MUST be larger or equal than N.
/// Use the set_* functions to set the data vector b and the matrix A.
/// Use the run_* functions in the correct order to solve your system (run_decomp_svd, then run_solve_svd)
/// You can score the result with run_score_svd_on_matrix
/// You can retrieve your solution with get_svd_solution.
///
/// @return
/// The score of the fitting : sqrt( ||A x - b||^2 ) with run_score_svd_on_matrix();
/// The fitted vector x with get_svd_solution();
///
/// @remarks
/// Calls in a wrong order of the functions will abort the program (i.e. if you try to solve the problem before you set a matrix A)
/// Once the matrix is decomposed, you can change the vector b and solve Ax=b with the new vector. (That's why those 2 functions are separated)
/// The matrix A is necessary to calculate the score (argument of run_score_svd_on_matrix), but the matrix A is not stored within
/// the SVD_solver object, so make sure you have it available when scoring (this is done on purpose for speed up)
/// Is it possible to speed up calculations by using FArraynD.index() call? ObjexxFCL doc is not really clear.
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


#ifndef INCLUDED_basic_svd_SVD_Solver_hh
#define INCLUDED_basic_svd_SVD_Solver_hh

// Package headers

#include <platform/types.hh>   // for Size
#include <string>              // for allocator
#include <utility/vector1.hh>  // for vector1

// Objexx headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

namespace basic {
namespace svd {

using namespace ObjexxFCL;

class SVD_Solver {

private:

	utility::vector1<double> cstyle_b_; //The vector b you want to use to solve Ax = b
	utility::vector1< utility::vector1<double> > cstyle_A_decomp_;
	utility::vector1< utility::vector1<double> > cstyle_v_;
	utility::vector1<double> cstyle_x_; //The solution x
	utility::vector1<double> cstyle_w_;
	utility::vector1<double> cstyle_tmp_;

	platform::Size M_; //M_ is the size of vector b (data)
	platform::Size N_; //N_ is the size of the parameter vector to optimize

	//A few states used to make sure the order of calls are correct
	bool b_is_set_;
	bool A_is_set_;
	bool A_is_decomp_;
	bool x_is_solved_;

public:

	SVD_Solver(); //construct

	~SVD_Solver(); //destruct

	SVD_Solver(platform::Size const M, platform::Size const N);

	SVD_Solver &
	operator=(SVD_Solver const & other); // =

	SVD_Solver(SVD_Solver const & other);

	/// @brief set the vector b of Ax=b
	void
	set_vector_b( utility::vector1<double>  const & b );

	/// @brief set the matrix A of Ax=b
	void
	set_matrix_A( utility::vector1< utility::vector1<double> >  const & A);

	/// @brief set the vector b of Ax=b (FArray1D version)
	void
	set_vector_b( FArray1D< double > const & b );

	/// @brief set the matrix A of Ax=b (FArray2D version)
	void
	set_matrix_A( FArray2D< double > const & A );


	/// @brief decompose the matrix A.
	/// Can be called after the matrix A and vector b are set (with set_matrix_A and set_vector_b)
	void
	run_decomp_svd();

	/// @brief minimize the cost sqrt( ||A x - b||^2 ) (but doesn't calculate it)
	/// Can be called after run_decomp_svd()
	void
	run_solve_svd();

	/// @brief return the score given the matrix A
	/// Can be called after run_decomp_svd()
	double
	run_score_svd_on_matrix(utility::vector1< utility::vector1<double> > const & cppstyle_A) const;

	/// @brief return the score given the matrix A (FArray2D version)
	/// Can be called after run_decomp_svd()
	double
	run_score_svd_on_matrix(FArray2D< double > const & A) const;

	/// @brief return the minimzed score without the need to call run_solve_svd()
	/// Can be called after run_decomp_svd()
	double
	run_score_svd_without_solving();

	/// @brief return the vector x that minimize ||Ax - b||^2
	/// Can be called after run_solve_svd()
	utility::vector1< double > const &
	get_svd_solution() const;

private:

	double
	pythag(double const & a, double const & b) const;

	void
	svdcmp();

	void
	svbksb();

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}//namespace svd
}//namespace basic


#endif
