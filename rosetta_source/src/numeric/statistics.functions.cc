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

} // namespace statistics
} // namespace numeric
