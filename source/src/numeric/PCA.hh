// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PCA.cc
///
/// @brief
/// @author Tim Jacobs
/// @author N-dimensional principal component analysis code by Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_numeric_PCA_hh
#define INCLUDED_numeric_PCA_hh

// Unit headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/types.hh>
//#include <numeric/linear_algebra/GeneralizedEigenSolver.hh>

// External headers
#include <Eigen/Dense>
#include <Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h>

// Utility headers
#include <utility/vector1.hh>

namespace numeric {

/// @brief return the first principal component
/// of the given set of points
template< typename T >
inline
xyzVector< T >
first_principal_component( utility::vector1< xyzVector< T > > const & coords ){
	return principal_components(coords).col(1);
}

/// @brief return a matrix containing the first 3 principal components
/// of the given set of points. Matrix columns are principal components,
/// first column is first component, etc.
template< typename T >
inline
xyzMatrix< T >
principal_components( utility::vector1< xyzVector< T > > const & coords ){
	return principal_components_and_eigenvalues(coords).first;
}

/// @brief return a vector containing the eigenvalues corresponding to the
/// first 3 principal components of the given set of points.
template< typename T >
inline
xyzVector< T >
principal_component_eigenvalues( utility::vector1< xyzVector< T > > const & coords ){
	return principal_components_and_eigenvalues(coords).second;
}

/// @brief return a pair containing a matrix of the first 3 principal components
/// and a vector of the corresponding eigenvalues of the given set of points.
template< typename T >
inline
std::pair<xyzMatrix<T>, xyzVector<T> >
principal_components_and_eigenvalues( utility::vector1< xyzVector< T > > const & coords )
{
	Size n_coords = coords.size();

	xyzVector< T > mean_vector(0.0);
	for ( Size i = 1; i <= n_coords; ++i ) {
		mean_vector += coords[i];
	}
	mean_vector /= n_coords;

	//Compute the covariance matrix
	xyzMatrix< T > covariance_matrix(0.0);
	for ( Size i = 1; i <= 3; ++i ) {
		for ( Size j = 1; j <= 3; ++j ) {
			for ( Size k = 1; k <= n_coords; ++k ) {
				covariance_matrix(i,j) += (mean_vector(i) - coords[k](i)) *
					(mean_vector(j) - coords[k](j));
			}
			covariance_matrix(i,j) /= n_coords;
		}
	}

	//Solve eigenvectors/values
	xyzMatrix< T > evecs;
	xyzVector< T > evals =
		eigenvector_jacobi( covariance_matrix, (T)0.001, evecs );

	utility::vector1< std::pair<T, xyzVector<T> > > sorted_val_vec_pairs;
	for ( Size i = 1 ; i <= 3; ++i ) {
		sorted_val_vec_pairs.push_back(std::make_pair(evals(i), evecs.col(i)));
	}

	//Sort eigenvectors highest to lowest based on eigenvalues
	sort(sorted_val_vec_pairs.rbegin(), sorted_val_vec_pairs.rend());
	xyzMatrix< T > sorted_evecs;
	sorted_evecs.col_x(sorted_val_vec_pairs[1].second);
	sorted_evecs.col_y(sorted_val_vec_pairs[2].second);
	sorted_evecs.col_z(sorted_val_vec_pairs[3].second);

	xyzVector< T > sorted_evals;
	sorted_evals.x(sorted_val_vec_pairs[1].first);
	sorted_evals.y(sorted_val_vec_pairs[2].first);
	sorted_evals.z(sorted_val_vec_pairs[3].first);

	std::pair<xyzMatrix<T>, xyzVector<T> > evecs_and_evals;
	evecs_and_evals = (std::make_pair(sorted_evecs, sorted_evals));

	return evecs_and_evals;
}

/// @brief Return a pair containing a matrix (vector of vectors) of all of the
/// principal components and a vector of the corresponding eigenvalues of the
/// given set of points in n-dimensional space.
/// @details Note that this does not assume that the input vectors are 3-dimensional.
/// If shift_center=false, the mean vector is not subtracted by this function.
/// (Failure to subtract mean vector prior to function call will produce odd results,
/// however.)
/// @author Vikram K. Mulligan (vmullig@uw.edu)
inline
std::pair<utility::vector1< utility::vector1< Real > >, utility::vector1< Real > >
principal_components_and_eigenvalues_ndimensions( utility::vector1< utility::vector1< Real > > const & coords, bool const shift_center )
{
	Size const n_coords = coords.size();
	runtime_assert_string_msg( n_coords > 0, "Empty coords matrix was passed to numeric::principal_components_and_eigenvalues_ndimensions()." );
	Size const n_dimensions = coords[1].size();

	utility::vector1< Real > mean_vector;
	mean_vector.resize(n_dimensions, 0);
	for ( Size i = 1; i <= n_coords; ++i ) {
		runtime_assert_string_msg( coords[i].size()==n_dimensions, "Not all coordinate vectors in the coords matrix passed to numeric::principal_components_and_eigenvalues_ndimensions() are of the same dimension." );
		if ( shift_center ) {
			for ( Size j=1; j<=n_coords; ++j ) {
				mean_vector[j] += coords[i][j];
			}
		}
	}
	if ( shift_center ) {
		for ( Size i=1; i<=n_dimensions; ++i ) {
			mean_vector[i] /= n_coords;
		}
	}

	//Compute the covariance matrix (upper-right half first):
	Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> covariance_matrix;
	covariance_matrix.resize(n_dimensions, n_dimensions);
	for ( Size i = 1; i <= n_dimensions; ++i ) {
		for ( Size j = i; j <= n_dimensions; ++j ) {
			covariance_matrix(i-1,j-1) = 0;
			for ( Size k = 1; k <= n_coords; ++k ) {
				covariance_matrix(i-1,j-1) += (mean_vector[i] - coords[k][i]) * (mean_vector[j] - coords[k][j]);
			}
			covariance_matrix(i-1,j-1) /= n_coords;
		}
	}
	//The covariance matrix is symmetric (so do the lower-left half by copying the upper-right half):
	for ( Size i=1; i<=n_dimensions; ++i ) {
		for ( Size j=1; j<i; ++j ) covariance_matrix(i-1,j-1)=covariance_matrix(j-1,i-1);
	}

	//Solve eigenvectors/values
	Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> >  eigensolver( covariance_matrix );
	Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> evecs = eigensolver.eigenvectors();
	Eigen::Matrix<Real, Eigen::Dynamic, 1> evals = eigensolver.eigenvalues();

	utility::vector1< std::pair<Real, utility::vector1<Real> > > sorted_val_vec_pairs;
	for ( Size i = 1 ; i <= n_dimensions; ++i ) {
		utility::vector1< Real > tempvect;
		tempvect.resize(n_dimensions);
		for ( Size j=1; j<=n_dimensions; ++j ) tempvect[j] = evecs(j-1,i-1);
		sorted_val_vec_pairs.push_back(std::make_pair(evals(i-1,0), tempvect));
	}

	//Sort eigenvectors highest to lowest based on eigenvalues
	sort(sorted_val_vec_pairs.rbegin(), sorted_val_vec_pairs.rend());
	utility::vector1< utility::vector1< Real > > sorted_evecs;
	sorted_evecs.resize(n_dimensions);
	for ( Size i=1; i<=n_dimensions; ++i ) {
		sorted_evecs[i]=sorted_val_vec_pairs[i].second;
	}

	utility::vector1 < Real > sorted_evals;
	sorted_evals.resize(n_dimensions);
	for ( Size i=1; i<=n_dimensions; ++i ) {
		sorted_evals[i] = sorted_val_vec_pairs[i].first;
	}

	std::pair<utility::vector1< utility::vector1< Real > >, utility::vector1< Real > > evecs_and_evals;
	evecs_and_evals = (std::make_pair(sorted_evecs, sorted_evals));

	return evecs_and_evals;
}

}//namespace

#endif  // INCLUDED_numeric_PCA_hh
