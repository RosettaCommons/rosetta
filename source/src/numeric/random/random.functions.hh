// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random.functions.hh
/// @brief  xyzVector and xyzMatrix functions
/// @author Justin R. Porter

#ifndef INCLUDED_numeric_random_random_functions_hh
#define INCLUDED_numeric_random_random_functions_hh

#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

namespace numeric {
namespace random {

/// @detail generate uniformly distributed vector on the unit sphere
/// @author Zhe Zhang
/// @author Justin R. Porter
// template< typename T >
// inline
// xyzVector< T >
// random_point_on_unit_sphere( RandomGenerator& RG ){
//   T phi = RG.uniform() * NumericTraits< T >::pi_2();
//   T theta = std::acos( sin_cos_range( 1-2*RG.uniform() ) );
//   return spherical_to_xyz( sphericalVector<T>( phi, theta, 1 ) );
// }

/// @detail generate uniformly distributed unit vector
/// @detail a random normal distribution of coordinates gives a uniform distribution of directions.
/// @detail same functionality as the commented random_point_on_unit_sphere(), but simpler algorithm
/// @detail rename as random_uniform_distributed_direction()?
/// @author Zhe Zhang
template< typename T >
inline
xyzVector< T >
random_point_on_unit_sphere( RandomGenerator& RG ){
	xyzVector< T > vec( RG.gaussian(), RG.gaussian(), RG.gaussian() );
	while ( vec.length() == 0 ) {
		vec = xyzVector< T >( RG.gaussian(), RG.gaussian(), RG.gaussian() );
	}
	return vec.normalized();
}

/// @detail generate axis and angle for axis-angle rotation for random rotation move in R/RT degrees of
/// freedom. rotation axis: uniformly distributed on unit sphere, rotation angle: chosen to mimic the
/// distribution of rotation angles obtained from gaussian distrbuted Euler angles (core/kinematics/Jump.cc),
/// which is a gamma-distribution-like distribution
/// @note gaussian distributed Euler angles do not give unbiased sampling in rotational space
/// by applying this angle to a uniformly chosen rotation axis unbiased rotational sampling is achieved
/// @brief gamma-distribution-like random angle generation, rot_mag makes exactly the same sense as in gaussian_move
/// @author Zhe Zhang
/// @auhtor Justin R. Porter

template < typename T >
inline
T // angle in radians
random_rotation_angle(
	T rot_mag,
	RandomGenerator & RG
)
{
	xyzMatrix< T > const mat(
		numeric::z_rotation_matrix_degrees( rot_mag * RG.gaussian() ) *
		numeric::y_rotation_matrix_degrees( rot_mag * RG.gaussian() ) *
		numeric::x_rotation_matrix_degrees( rot_mag * RG.gaussian() ) );
	T theta;
	rotation_axis( mat, theta );
	return theta;
}

/// @brief produce a random translation in xyz space
/// @author Zhe Zhang
/// @author Justin R. Porter
template < typename T >
inline
xyzVector< T >
random_translation(
	T trans_mag,
	RandomGenerator& RG
)
{
	return xyzVector<T>(
		trans_mag * RG.gaussian(),
		trans_mag * RG.gaussian(),
		trans_mag * RG.gaussian() );
}

/// @brief Return the index in the CDF array such that it is smaller than or equal to the
/// uniform-random-number (urn) and that the next entry in the array is larger than urn.
/// The CDF ought to represent the exclusive cumulative sum of the probabilities of some
/// discrete set so that the width for entry i -- the gap between entry i and entry i+1
/// should be the probability of entry i.
template < typename T >
platform::Size
binary_search_cdf(
	utility::vector1< T > const & cdf,
	T urn
)
{
	platform::Size lower( 1 ), upper( cdf.size() );
	platform::Size guess( 0 );

	while ( true ) {

		guess = ( upper + lower ) / 2;
		T guess_val = cdf[ guess ];
		T next_val = guess < cdf.size() ? cdf[ guess+1 ] : 1.1;
		if ( guess_val <= urn && urn < next_val ) {
			// found it!
			break;
		}
		if ( guess_val <= urn ) {
			lower = guess+1;
		} else {
			upper = guess-1;
		}
	}

	return guess;

}

/// @brief Choose an element from an array of positions representing the
/// cumulative distribution function for some target distribution.  This uses
/// binary search to quickly find the target index.  If any position has the
/// same probability as the preceeding position, then its index won't be
/// returned.
template < typename T >
platform::Size
pick_random_index_from_cdf(
	utility::vector1< T > const & cdf,
	RandomGenerator & RG
)
{
	T const urn = RG.uniform();
	return binary_search_cdf( cdf, urn );
}


}
}

#endif
