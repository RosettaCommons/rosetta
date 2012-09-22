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
#include <numeric/types.hh>

// External headers
#include <Eigen/Eigen>

namespace numeric {

	/// @brief a function to sort pairs that doesn't require
	/// operators for the second element in the pair. This is
	/// needed below
	inline bool
	compare_first_only(const std::pair<double,Eigen::VectorXd> &left, const std::pair<double,Eigen::VectorXd> &right) {
			return left.first < right.first;
	}

	/// @brief return the first principal component
	/// of the given set of points
	template< typename T >
	inline
	xyzVector< T >
	first_principal_component( utility::vector1< xyzVector< T > > const & coords ){
		return principal_components(coords).col(1);
	}

	/// @brief return a matrix containing the 3 principal components
	/// of the given set of points
	template< typename T >
	inline
	xyzMatrix< T >
	principal_components( utility::vector1< xyzVector< T > > const & coords )
	{
		using namespace Eigen;

		Size n = coords.size();
		MatrixXd data_points = MatrixXd::Zero(3, n);//3 dimensions x # of points
		for(Size i=1; i<=coords.size(); ++i)
		{
			data_points(0,i-1)=coords[i].x();
			data_points(1,i-1)=coords[i].y();
			data_points(2,i-1)=coords[i].z();
		}

		MatrixXd mean_subtracted_data = data_points;
		for (int i = 0; i < 3; ++i)
		{
			T mean = (mean_subtracted_data.row(i).sum())/n; //compute mean of each dimension (x,y,z)
			VectorXd meanVector  = VectorXd::Constant(n,mean); //create a vector with constant value = mean
			mean_subtracted_data.row(i) -= meanVector; //subtract mean from every point for the current dimension
		}

		// get the covariance matrix
		MatrixXd Covariance = MatrixXd::Zero(3, 3);
		Covariance = (1 / (T)n) * mean_subtracted_data * mean_subtracted_data.transpose();

		// compute the eigenvalue on the Cov Matrix
		EigenSolver<MatrixXd> m_solve(Covariance);
		VectorXd eigenvalues = VectorXd::Zero(3);
		eigenvalues = m_solve.eigenvalues().real();

		MatrixXd eigenvectors = MatrixXd::Zero(n, 3);
		eigenvectors = m_solve.eigenvectors().real();

		//Create a mapping between eigenvalues and eigenvectors
		utility::vector1< std::pair<double, VectorXd > > value_vector_pairs;
		for (Size i = 1 ; i <= 3; ++i)
		{
			value_vector_pairs.push_back(std::make_pair(eigenvalues(i-1), eigenvectors.col(i-1)));//eigenvalues is 0-indexed
		}

		//Sort eigenvectors highest to lowest based on eigenvalues
		sort(value_vector_pairs.begin(), value_vector_pairs.end(), compare_first_only);

		//Convert to xyzMatrix
		xyzMatrix<T> sorted_eigenvectors;
		for (int i = 0; i < 3; i++)
		{
			numeric::xyzVector<T> xyz_eigenvector(
				value_vector_pairs[i+1].second(0),
				value_vector_pairs[i+1].second(1),
				value_vector_pairs[i+1].second(2));
			sorted_eigenvectors.col(i+1, xyz_eigenvector);
		}

		return sorted_eigenvectors;
	}

}//namespace

#endif  // INCLUDED_numeric_PCA_hh
