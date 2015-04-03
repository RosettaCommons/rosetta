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
 /// @file basic/svd/SVD_Solver.cc
 ///
 /// @brief SVD solver class
 ///`
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
 /// The fited vector x with get_svd_solution();
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
 /// @authorv Christophe Schmitz & Srivatsan Raman
 ///
 ////////////////////////////////////////////////


// Unit headers
#include <basic/svd/SVD_Solver.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Basic headers

// Objexx headers
#include <ObjexxFCL/Fmath.hh>

// C++ headers
#include <iostream>

#include <platform/types.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

namespace basic{
namespace svd{

SVD_Solver::~SVD_Solver(){
}

SVD_Solver &
SVD_Solver::operator=(SVD_Solver const & other){

	/*
	if ((M_ != other.M_) || (N_ != other.N_)){
		//utility_exit_with_message( "You can't call the = operator on SVD_Solver object of different size" );

		M_ = other.M_;
		N_ = other.N_;

		if(N_ >= M_){
			utility_exit_with_message("First parameter of SVD_Solver constructor MUST be larger than the second parameter");
		}

		cstyle_A_decomp_.resize(other.M_, other.N_);
		cstyle_b_.resize(other.M_);
		cstyle_v_.resize(other.N_, other.N_);
		cstyle_x_.resize(other.N_);
		cstyle_w_.resize(other.N_);
		cstyle_tmp_.resize(other.N_);

		cstyle_b_ = other.cstyle_b_;
		cstyle_A_decomp_ = other.cstyle_A_decomp_;
		cstyle_v_ = other.cstyle_v_;
		cstyle_x_ = other.cstyle_x_;
		cstyle_w_ = other.cstyle_w_;
		cstyle_tmp_ = other.cstyle_tmp_;

		b_is_set_ = other.b_is_set_;
		A_is_set_ = other.A_is_set_;
		A_is_decomp_ = other.A_is_decomp_;
		x_is_solved_ = other.x_is_solved_;


	} else if ( this != &other ) {
		cstyle_b_ = other.cstyle_b_;
		cstyle_A_decomp_ = other.cstyle_A_decomp_;
		cstyle_v_ = other.cstyle_v_;
		cstyle_x_ = other.cstyle_x_;
		cstyle_w_ = other.cstyle_w_;
		cstyle_tmp_ = other.cstyle_tmp_;

		b_is_set_ = other.b_is_set_;
		A_is_set_ = other.A_is_set_;
		A_is_decomp_ = other.A_is_decomp_;
		x_is_solved_ = other.x_is_solved_;
	}
	return *this;
	*/

	if ( this != &other ) {
		M_ = other.M_;
		N_ = other.N_;

		if(N_ >= M_){
			utility_exit_with_message("First parameter of SVD_Solver constructor MUST be larger than the second parameter");
		}

		platform::Size i;
		cstyle_A_decomp_.resize(other.M_);
		for (i = 1; i <= other.M_; ++i){
			cstyle_A_decomp_[i].resize(other.N_);
		}
		cstyle_b_.resize(other.M_);

		cstyle_v_.resize(other.N_);
		for (i = 1; i <= other.N_; ++i){
			cstyle_v_[i].resize(other.N_);
		}
		cstyle_x_.resize(other.N_);
		cstyle_w_.resize(other.N_);
		cstyle_tmp_.resize(other.N_);

		cstyle_b_ = other.cstyle_b_;
		cstyle_A_decomp_ = other.cstyle_A_decomp_;
		cstyle_v_ = other.cstyle_v_;
		cstyle_x_ = other.cstyle_x_;
		cstyle_w_ = other.cstyle_w_;
		cstyle_tmp_ = other.cstyle_tmp_;

		b_is_set_ = other.b_is_set_;
		A_is_set_ = other.A_is_set_;
		A_is_decomp_ = other.A_is_decomp_;
		x_is_solved_ = other.x_is_solved_;
	}

	return *this;
}

SVD_Solver:: SVD_Solver(SVD_Solver const & other)//:
	//	M_(other.M_), N_(other.N_)
{

	M_ = other.M_;
	N_ = other.N_;

  cstyle_b_ = other.cstyle_b_;
  cstyle_A_decomp_ = other.cstyle_A_decomp_;
  cstyle_v_ = other.cstyle_v_;
  cstyle_x_ = other.cstyle_x_;
  cstyle_w_ = other.cstyle_w_;
  cstyle_tmp_ = other.cstyle_tmp_;

	b_is_set_ = other.b_is_set_;
	A_is_set_ = other.A_is_set_;
	A_is_decomp_ = other.A_is_decomp_;
	x_is_solved_ = other.x_is_solved_;
}

SVD_Solver::SVD_Solver()//:
	//	M_(0), N_(0)
{
	utility_exit_with_message( "You shouldn't call the empty constructor for SVD_Solver class" );
}

///////////////////////////////////////////////
/// @brief M is the size of vector b (experimental data); N is the size of the vector to optimize M >= N
///////////////////////////////////////////////
SVD_Solver::SVD_Solver(platform::Size const M, platform::Size const N)//:
//	M_(M), N_(N)
{

	M_ = M;
	N_ = N;

	if(N_ >= M_){
		utility_exit_with_message("First parameter of SVD_Solver constructor MUST be larger than the second parameter");
	}
	platform::Size i;

	cstyle_A_decomp_.resize(M);
	for (i = 1; i <= M; ++i){
		cstyle_A_decomp_[i].resize(N);
	}

	cstyle_b_.resize(M);

	cstyle_v_.resize(N);
	for (i = 1; i <= N; ++i){
		cstyle_v_[i].resize(N);
	}

	cstyle_x_.resize(N);
	cstyle_w_.resize(N);
	cstyle_tmp_.resize(N);

	b_is_set_ = false;
	A_is_set_ = false;
	A_is_decomp_ = false;
	x_is_solved_ = false;
}

///////////////////////////////////////////////
/// @brief To minimize ||Ax - b||^2, you need to set the vector b
///////////////////////////////////////////////
void
SVD_Solver::set_vector_b(utility::vector1<double>  const & b){
	platform::Size i;
	if (M_ != b.size()){
		utility_exit_with_message("Size dimension differs when trying to set vector b in SVD_solver class");
	}

	for (i = 1; i <= M_; ++i){
		cstyle_b_[i] = b[i];
	}

	b_is_set_ = true;
	x_is_solved_ = false;
}

///////////////////////////////////////////////
/// @brief To minimize ||Ax - b||^2, you need to set the vector b
///////////////////////////////////////////////
void
SVD_Solver::set_vector_b(FArray1D< double > const & b){
	platform::Size i;
	if (M_ != b.size()){
		utility_exit_with_message("Size dimension differs when trying to set vector b in SVD_solver class");
	}

	for (i = 1; i <= M_; ++i){
		cstyle_b_[i] = b(i);
	}

	b_is_set_ = true;
	x_is_solved_ = false;
}


///////////////////////////////////////////////
/// @brief To minimize ||Ax - b||^2, you need to set the matrix A
///////////////////////////////////////////////
void
SVD_Solver::set_matrix_A(utility::vector1< utility::vector1<double> >  const & A){
	platform::Size i, j;
	if (M_ != A.size()){
		utility_exit_with_message("Size dimension differs when trying to set the matrix A in SVD_solver class");
	}

	for (i = 1; i <= M_; ++i){
		if (N_ != A[i].size()){
			utility_exit_with_message( "Size dimension differs when trying to set the matrix A in SVD_solver class");
		}
		for (j = 1; j <= N_; ++j){
			cstyle_A_decomp_[i][j] = A[i][j];
		}
	}

	A_is_set_ = true;
	A_is_decomp_ = false;
	x_is_solved_ = false;
}


///////////////////////////////////////////////
/// @brief To minimize ||Ax - b||^2, you need to set the matrix A
///////////////////////////////////////////////
void
SVD_Solver::set_matrix_A( FArray2D< double > const & A){
	platform::Size i, j;

	if (M_*N_ != A.size()){
		utility_exit_with_message( "Size dimension differs when trying to set the matrix A in SVD_solver class");
	}

	for (i = 1; i <= M_; ++i){
		for (j = 1; j <= N_; ++j){
			cstyle_A_decomp_[i][j] = A(i, j);
		}
	}

	A_is_set_ = true;
	A_is_decomp_ = false;
	x_is_solved_ = false;
}

///////////////////////////////////////////////
/// @brief Decomposition of the matrix A. Can be called only when A is set.
/// You can't decompose twice in a row without reseting the matrix A.
///////////////////////////////////////////////
void
SVD_Solver::run_decomp_svd(){
	if(!A_is_set_){
		utility_exit_with_message("The matrix A is not ready to be decomposed in SVD_solver class");
	}
	svdcmp();
	A_is_decomp_ = true;
}

///////////////////////////////////////////////
/// @brief Minimize ||Ax - b||^2. Can be called only when b is set and A is decomposed.
/// You are allowed to update the vector b and minimize again.
///////////////////////////////////////////////
void
SVD_Solver::run_solve_svd(){
	if((!A_is_decomp_) || (!b_is_set_)){
		utility_exit_with_message("SVD_solver object is not in a state to solve Ax = b");
	}
	svbksb();
	x_is_solved_ = true;

}

///////////////////////////////////////////////
/// @brief Return the optimized vector x. Can be called only when ||Ax-b||^2 has been minimized.
///////////////////////////////////////////////
/*
FArray1D< double > const &
SVD_Solver::get_svd_solution() const{
	if(!x_is_solved_){
		utility_exit_with_message("SVD_solver object has not yet solved Ax = b");
	}
	return(cstyle_x_);
}
*/

utility::vector1< double > const &
SVD_Solver::get_svd_solution() const{
	if(!x_is_solved_){
		utility_exit_with_message("SVD_solver object has not yet solved Ax = b");
	}
	return(cstyle_x_);
}


///////////////////////////////////////////////
/// @brief Return the score SQRT(||Ax-b||^2). Can be called only when ||Ax-b||^2 has been minimized.
/// You need the give the original A matrix as argument.
///////////////////////////////////////////////
double
SVD_Solver::run_score_svd_on_matrix(utility::vector1< utility::vector1<double> > const & cppstyle_A) const{

	platform::Size i, j;
	double score, temp;

	if(!x_is_solved_){
		utility_exit_with_message("SVD_solver object is not in a state to score ||Ax = b||^2");
	}
	score = 0;
	for(i = 1; i <= M_; ++i){
		temp = 0;
		for(j = 1; j <= N_; ++j){
			temp += cppstyle_A[i][j] * cstyle_x_[j];
		}
		score += ( temp - cstyle_b_[i] ) * ( temp - cstyle_b_[i] ) ;
	}


	return(sqrt(score));
}

///////////////////////////////////////////////
/// @brief Return the score SQRT(||Ax-b||^2). Can be called only when ||Ax-b||^2 has been minimized.
/// You need the give the original A matrix as argument.
///////////////////////////////////////////////
double
SVD_Solver::run_score_svd_on_matrix(FArray2D< double > const & A) const{

	platform::Size i, j;
	double score, temp;

	if(!x_is_solved_){
		utility_exit_with_message("SVD_solver object is not in a state to score ||Ax = b||^2");
	}
	score = 0;
	for(i = 1; i <= M_; ++i){
		temp = 0;
		for(j = 1; j <= N_; ++j){
			temp += A(i, j) * cstyle_x_[j];
		}
		score += ( temp - cstyle_b_[i] ) * ( temp - cstyle_b_[i] ) ;
	}


	return(sqrt(score));
}

///////////////////////////////////////////////
/// @brief Attempt to speed up calculation of the cost without actually solving Ax = b
/// The routine works, but is not faster (nor slower)
///////////////////////////////////////////////
double
SVD_Solver::run_score_svd_without_solving(){

	platform::Size i, j;
	double score, temp;

	if((!A_is_decomp_) || (!b_is_set_)){
		utility_exit_with_message("SVD_Solver object not in state to call run_score_svd_without_solving");
	}

	for(i = 1; i <= N_; ++i){
		cstyle_tmp_[i] = 0;
		for(j = 1; j <= M_; j++){
			cstyle_tmp_[i] += cstyle_A_decomp_[j][i] *  cstyle_b_[j];
		}
	}

	score = 0;
	for(i = 1; i <= M_; ++i){
		temp = 0;
		for(j = 1; j <= N_; ++j){
			temp += cstyle_A_decomp_[i][j] * cstyle_tmp_[j];
		}
		temp -= cstyle_b_[i];
		score += 	temp * temp;
	}

	return(sqrt(score));
}

////////////////////////////////////////////////////////////////
// copy past from ResidualDipolarCouplingEnergy.cc of Srivatsan Raman
////////////////////////////////////////////////////////////////
double
SVD_Solver::pythag(double const & a, double const & b) const {

	double pythag;

	double absa = std::abs(a);
	double absb = std::abs(b);
	if ( absa > absb ) {
		double const ratio = absb/absa;
		pythag = absa * std::sqrt( 1.0 + ( ratio * ratio ) );
	} else {
		if ( absb == 0.0 ) {
			pythag = 0.0;
		} else {
			double const ratio = absa/absb;
			pythag = absb * std::sqrt( 1.0 + ( ratio * ratio ) );
		}
	}
	return pythag;
}


////////////////////////////////////////////////////////////////
// Copy past from ResidualDipolarCouplingEnergy.cc of Srivatsan Raman
// Can this be optimized with FArraynD.index() call?
// Is it somehhow possible to score quicker without calculating the x vector?
////////////////////////////////////////////////////////////////
void
SVD_Solver::svbksb(){
	double s;

	for ( platform::Size j = 1; j <= N_; ++j ) {
		s = 0.0;
		if ( cstyle_w_[j] != 0.0 ) {
			for ( platform::Size i = 1; i <= M_; ++i ) {
				s += cstyle_A_decomp_[i][j] * cstyle_b_[i];
			}
			s /= cstyle_w_[j];
		}
		cstyle_tmp_[j] = s;
	}
	for ( platform::Size j = 1; j <= N_; ++j ) {
		s = 0.0;
		for ( platform::Size jj = 1; jj <= N_; ++jj ) {
			s += cstyle_v_[j][jj] * cstyle_tmp_[jj];
		}
		cstyle_x_[j] = s;
	}
}

////////////////////////////////////////////////////////////////
// copy past from ResidualDipolarCouplingEnergy.cc of Srivatsan Raman
// Can this be optimized with FArraynD.index() call?
////////////////////////////////////////////////////////////////
void
SVD_Solver::svdcmp(){

//U    USES pythag
	platform::Size i,its,j,jj,k,l,nm;
	double anorm, c, f, g, h, s, scale, x, y, z;
	g = 0.0;
	scale = 0.0;
	anorm = 0.0;
	l = 0;
	nm = 0;
	for ( i = 1; i <= N_; ++i ) {
		l = i+1;
		cstyle_tmp_[i] = scale*g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( i <= M_ ) {
			for ( k = i; k <= M_; ++k ) {
				scale += std::abs(cstyle_A_decomp_[k][i]);
			}
			if ( scale != 0.0 ) {
				for ( k = i; k <= M_; ++k ) {
					cstyle_A_decomp_[k][i] /= scale;
					s += cstyle_A_decomp_[k][i]*cstyle_A_decomp_[k][i];
				}
				f = cstyle_A_decomp_[i][i];
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				cstyle_A_decomp_[i][i] = f-g;
				for ( j = l; j <= N_; ++j ) {
					s = 0.0;
					for ( k = i; k <= M_; ++k ) {
						s += cstyle_A_decomp_[k][i]*cstyle_A_decomp_[k][j];
					}
					f = s/h;
					for ( k = i; k <= M_; ++k ) {
						cstyle_A_decomp_[k][j] += f*cstyle_A_decomp_[k][i];
					}
				}
				for ( k = i; k <= M_; ++k ) {
					cstyle_A_decomp_[k][i] *= scale;
				}
			}
		}
		cstyle_w_[i] = scale *g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( (i <= M_) && (i != N_) ) {
			for ( k = l; k <= N_; ++k ) {
				scale += std::abs(cstyle_A_decomp_[i][k]);
			}
			if ( scale != 0.0 ) {
				for ( k = l; k <= N_; ++k ) {
					cstyle_A_decomp_[i][k] /= scale;
					s += cstyle_A_decomp_[i][k]*cstyle_A_decomp_[i][k];
				}
				f = cstyle_A_decomp_[i][l];
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				cstyle_A_decomp_[i][l] = f-g;
				for ( k = l; k <= N_; ++k ) {
					cstyle_tmp_[k] = cstyle_A_decomp_[i][k]/h;
				}
				for ( j = l; j <= M_; ++j ) {
					s = 0.0;
					for ( k = l; k <= N_; ++k ) {
						s += cstyle_A_decomp_[j][k]*cstyle_A_decomp_[i][k];
					}
					for ( k = l; k <= N_; ++k ) {
						cstyle_A_decomp_[j][k] += s*cstyle_tmp_[k];
					}
				}
				for ( k = l; k <= N_; ++k ) {
					cstyle_A_decomp_[i][k] *= scale;
				}
			}
		}
		anorm = std::max(anorm,(std::abs(cstyle_w_[i])+std::abs(cstyle_tmp_[i])));
	}
	for ( i = N_; i >= 1; --i ) {
		if ( i < N_ ) {
			if ( g != 0.0 ) {
				for ( j = l; j <= N_; ++j ) {
					cstyle_v_[j][i] = (cstyle_A_decomp_[i][j]/cstyle_A_decomp_[i][l])/g;
				}
				for ( j = l; j <= N_; ++j ) {
					s = 0.0;
					for ( k = l; k <= N_; ++k ) {
						s += cstyle_A_decomp_[i][k]*cstyle_v_[k][j];
					}
					for ( k = l; k <= N_; ++k ) {
						cstyle_v_[k][j] += s*cstyle_v_[k][i];
					}
				}
			}
			for ( j = l; j <= N_; ++j ) {
				cstyle_v_[i][j] = 0.0;
				cstyle_v_[j][i] = 0.0;
			}
		}
		cstyle_v_[i][i] = 1.0;
		g = cstyle_tmp_[i];
		l = i;
	}
	for ( i = std::min(M_,N_); i >= 1; --i ) {
		l = i+1;
		g = cstyle_w_[i];
		for ( j = l; j <= N_; ++j ) {
			cstyle_A_decomp_[i][j] = 0.0;
		}
		if ( g != 0.0 ) {
			g = 1.0/g;
			for ( j = l; j <= N_; ++j ) {
				s = 0.0;
				for ( k = l; k <= M_; ++k ) {
					s += cstyle_A_decomp_[k][i]*cstyle_A_decomp_[k][j];
				}
				f = (s/cstyle_A_decomp_[i][i])*g;
				for ( k = i; k <= M_; ++k ) {
					cstyle_A_decomp_[k][j] += f*cstyle_A_decomp_[k][i];
				}
			}
			for ( j = i; j <= M_; ++j ) {
				cstyle_A_decomp_[j][i] *= g;
			}
		} else {
			for ( j = i; j <= M_; ++j ) {
				cstyle_A_decomp_[j][i] = 0.0;
			}
		}
		cstyle_A_decomp_[i][i] += 1.0;
	}
	for ( k = N_; k >= 1; --k ) {
		for ( its = 1; its <= 30; ++its ) {
			for ( l = k; l >= 1; --l ) {
				nm = l-1;
				if ( (std::abs(cstyle_tmp_[l])+anorm) == anorm ) goto L2;
				if ( (std::abs(cstyle_w_[nm])+anorm) == anorm ) break;
			}
			c = 0.0;
			s = 1.0;
			for ( i = l; i <= k; ++i ) {
				f = s*cstyle_tmp_[i];
				cstyle_tmp_[i] *= c;
				if ( (std::abs(f)+anorm) == anorm ) break;
				g = cstyle_w_[i];
				h = pythag(f, g);
				cstyle_w_[i] = h;
				h = 1.0/h;
				c = (g*h);
				s = -(f*h);
				for ( j = 1; j <= M_; ++j ) {
					y = cstyle_A_decomp_[j][nm];
					z = cstyle_A_decomp_[j][i];
					cstyle_A_decomp_[j][nm] = (y*c)+(z*s);
					cstyle_A_decomp_[j][i] = -(y*s)+(z*c);
				}
			}
L2:
			z = cstyle_w_[k];
			if ( l == k ) {
				if ( z < 0.0 ) {
					cstyle_w_[k] = -z;
					for ( j = 1; j <= N_; ++j ) {
						cstyle_v_[j][k] = -cstyle_v_[j][k];
					}
				}
				break;
			}
			if ( its == 30) utility_exit_with_message("no convergence in svdcmp \n" );
			x = cstyle_w_[l];
			nm = k-1;
			y = cstyle_w_[nm];
			g = cstyle_tmp_[nm];
			h = cstyle_tmp_[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c = 1.0;
			s = 1.0;
			for ( j = l; j <= nm; ++j ) {
				i = j+1;
				g = cstyle_tmp_[i];
				y = cstyle_w_[i];
				h = s*g;
				g *= c;
				z = pythag(f,h);
				cstyle_tmp_[j] = z;
				c = f/z;
				s = h/z;
				f = (x*c)+(g*s);
				g = -(x*s)+(g*c);
				h = y*s;
				y *= c;
				for ( jj = 1; jj <= N_; ++jj ) {
					x = cstyle_v_[jj][j];
					z = cstyle_v_[jj][i];
					cstyle_v_[jj][j] = (x*c)+(z*s);
					cstyle_v_[jj][i] = -(x*s)+(z*c);
				}
				z = pythag(f,h);
				cstyle_w_[j] = z;
				if ( z != 0.0 ) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = (c*g)+(s*y);
				x = -(s*g)+(c*y);
				for ( jj = 1; jj <= M_; ++jj ) {
					y = cstyle_A_decomp_[jj][j];
					z = cstyle_A_decomp_[jj][i];
					cstyle_A_decomp_[jj][j] = (y*c)+(z*s);
					cstyle_A_decomp_[jj][i] = -(y*s)+(z*c);
				}
			}
			cstyle_tmp_[l] = 0.0;
			cstyle_tmp_[k] = f;
			cstyle_w_[k] = x;
		}
	}
}

}//namespace svd
}//namespace basic
