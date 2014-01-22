// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
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

#ifndef INCLUDED_numeric_PCA_hh
#define INCLUDED_numeric_PCA_hh

// Unit headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/types.hh>

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
		for (Size i = 1; i <= n_coords; ++i)
		{
			mean_vector += coords[i];
		}
		mean_vector /= n_coords;
		
		//Compute the covariance matrix
		xyzMatrix< T > covariance_matrix(0.0);
		for (Size i = 1; i <= 3; ++i)
		{
			for (Size j = 1; j <= 3; ++j)
			{
				for (Size k = 1; k <= n_coords; ++k)
				{
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
		for (Size i = 1 ; i <= 3; ++i)
		{
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

	
}//namespace

#endif  // INCLUDED_numeric_PCA_hh
