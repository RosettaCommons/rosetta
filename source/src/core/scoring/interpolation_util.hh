// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/interpolation/interpolation_util.hh
/// @brief  Useful function for interpolation...
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_interpolation_util_hh
#define INCLUDED_core_scoring_interpolation_util_hh

#include <core/types.hh>
#include <ObjexxFCL/FArray1D.hh>

namespace core {
namespace scoring {

/////////////////////////////////////////////////////////////////////////////////////////
// A generally useful function for interpolation.
// Its a little puzzling that this isn't in any of the utility::interpolation
// files, but I sure can't find it.
inline
void
interpolate_value_and_deriv(
	ObjexxFCL::FArray1D <Real> const & potential,
	Real const & bin_width,
	Real const & r,
	Real & value,
	Real & deriv )
{

	Size const & NUM_BINS( potential.size() );

	Size r_bin_lo = Size( r/bin_width + 0.5 );
	Size r_bin_hi  = r_bin_lo + 1;

	value = 0.0;
	deriv = 0.0;

	if ( r_bin_lo > NUM_BINS ) return;

	if ( r_bin_hi > NUM_BINS ) {
		value = potential( NUM_BINS );
	} else if ( r_bin_lo < 1 ) {
		value = potential( 1 );
	} else {
		Real const fraction = (r/bin_width + 0.5) - r_bin_lo;
		value = ( (1 - fraction ) * potential( r_bin_lo ) +
			fraction * potential( r_bin_hi ) );
		deriv = ( potential( r_bin_hi ) - potential( r_bin_lo ) ) / bin_width;
	}

}

} // end namespace scoring
} // end namespace core

#endif
