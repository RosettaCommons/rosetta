// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_ScoringUtil.cc
/// @author Rhiju Das

// Unit headers
#include <core/scoring/rna/RNA_ScoringUtil.hh>
#include <core/types.hh>
#include <cassert>

namespace core {
namespace scoring {
namespace rna {
///////////////////////////////////////////////////////////////////////////////
// Simple cubic spline.
void
get_fade_correction(
   Real const z,
	 Real const cutoff_lower,
	 Real const cutoff_upper,
	 Real const fade_zone,
	 Real & fade_value,
	 Real & fade_deriv )
{
	assert( fade_zone > 0 );

	fade_value = 1.0;
	fade_deriv = 0.0;

	if ( z < cutoff_lower || z > cutoff_upper ){
		fade_value = 0.0;
	} else if ( z < cutoff_lower + fade_zone ) {
		//Check little strip near lower cutoff.
		Real const b = -1.0 * ( z - ( cutoff_lower + fade_zone ) )/ fade_zone;
		Real const b2 = b*b;
		Real const b3 = b2*b;
		fade_value = ( 2 * b3 - 3 * b2 + 1 );
		fade_deriv = -1.0 * ( 6 * b2 - 6 * b ) / fade_zone;
	} else if ( z > cutoff_upper - fade_zone ) {
		//Check little strip near upper cutoff.
		Real const b =  ( z - ( cutoff_upper - fade_zone ) )/ fade_zone;
		Real const b2 = b*b;
		Real const b3 = b2*b;
		fade_value = ( 2 * b3 - 3 * b2 + 1 );
		fade_deriv = ( 6 * b2 - 6 * b ) / fade_zone;
	}

	return;

}

} //rna
} //scoring
} //core


