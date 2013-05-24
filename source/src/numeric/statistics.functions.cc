// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/statistics.functions.cc
/// @brief
/// @author James Thompson

// Platform headers
#include <platform/types.hh>
#include <numeric/statistics.functions.hh>
#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <cmath>

namespace numeric {
namespace statistics {

numeric::Real kl_divergence(
	utility::vector1< numeric::Real > const & prior,
	utility::vector1< numeric::Real > const & posterior
) {
	assert( prior.size() == posterior.size() );

	// p = prior distribution, q = posterior distribution
	// relative entropy = sum( p(i) * log( p(i) / q(i) )
	Real div( 0.0 );
	typedef utility::vector1< numeric::Real >::const_iterator iter;
	for ( iter p = prior.begin(), p_end = prior.end(),
				q = posterior.begin(), q_end = posterior.end();
				p != p_end && q != q_end; ++p, ++q
	) {
		div += *p * std::log( *p / *q );
	}

	return div;
}

numeric::Real corrcoef(
				 utility::vector1< numeric::Real > const & vec1,
				 utility::vector1< numeric::Real > const & vec2)

{
	assert( vec1.size() == vec2.size() );
	numeric::Real m1 = mean( vec1.begin(),vec1.end(),0.0 );
	numeric::Real m2 = mean( vec2.begin(),vec2.end(),0.0 );
	numeric::Real sd1 = std_dev_with_provided_mean( vec1.begin(), vec1.end(), m1 );
	numeric::Real sd2 = std_dev_with_provided_mean( vec2.begin(), vec2.end(), m2 );
	return corrcoef_with_provided_mean_and_std_dev( vec1, m1, sd1, vec2, m2, sd2 );
 }

numeric::Real corrcoef_with_provided_mean_and_std_dev(utility::vector1< numeric::Real > const & vec1,
																											numeric::Real m1,
																											numeric::Real sd1,
																											utility::vector1< numeric::Real > const & vec2,
																											numeric::Real m2,
																											numeric::Real sd2) 
{
	numeric::Real cov(0);
	Size n(0);
	typedef utility::vector1< numeric::Real >::const_iterator iter;
	for(iter v1=vec1.begin(),v2=vec2.begin(),v1_end=vec1.end(),v2_end=vec2.end();
			v1 != v1_end && v2 != v2_end;
			++v1,++v2,n++)
		cov+=(*v1-m1)*(*v2-m2);
	cov/=(n-1);
	return cov/(sd1*sd2); 
}
numeric::Real cov(
				 utility::vector1< numeric::Real > const & vec1,
				 utility::vector1< numeric::Real > const & vec2)

{
	assert( vec1.size() == vec2.size() );
	numeric::Real m1 = mean( vec1.begin(),vec1.end(),0.0 );
	numeric::Real m2 = mean( vec2.begin(),vec2.end(),0.0 );
	return cov_with_provided_mean( vec1, m1, vec2, m2 );
 }

numeric::Real cov_with_provided_mean(utility::vector1< numeric::Real > const & vec1,
									numeric::Real m1,
									utility::vector1< numeric::Real > const & vec2,
									numeric::Real m2)
{
	numeric::Real cov(0);
	Size n(0);
	typedef utility::vector1< numeric::Real >::const_iterator iter;
	for(iter v1=vec1.begin(),v2=vec2.begin(),v1_end=vec1.end(),v2_end=vec2.end();
			v1 != v1_end && v2 != v2_end;
			++v1,++v2,n++)
		cov+=(*v1-m1)*(*v2-m2);
	cov/=(n-1);
	return cov;
}
} // namespace statistics
} // namespace numeric
