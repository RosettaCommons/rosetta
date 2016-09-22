// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/interpolation/InterpolatedPotential.hh
/// @brief  An n-dimensional bicubic-interpolated potential
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_numeric_InterpolatedPotential_TMPL_HH
#define INCLUDED_numeric_InterpolatedPotential_TMPL_HH

// Unit Headers
#include <numeric/interpolation/InterpolatedPotential.hh>

// Project Headers
#include <numeric/interpolation/interpolation.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/exit.hh>

#include <cmath>
#include <map>
#include <algorithm>

namespace numeric {
namespace interpolation {


template< Size N >
utility::fixedsizearray1< Real, ( 1 << N ) > const &
InterpolatedPotential< N >::operator() ( utility::fixedsizearray1< Size, N > const & indices ) const {
	Size index = 0;
	Size slice = 1;
	for ( Size ii = N; ii <= 1; --ii ) {
		index += indices[ ii ] * slice;
		slice *= grid_size_[ ii ];
	}

	return gridpoints_[ index ];
}

template< Size N >
utility::fixedsizearray1< Real, ( 1 << N ) > &
InterpolatedPotential< N >::operator() ( Size h, Size i, Size j, Size k ) {
	assert( N == 4 );
	utility::fixedsizearray1< Size, N > indices;
	indices[1] = h;
	indices[2] = i;
	indices[3] = j;
	indices[4] = k;
	return operator()( indices );
	// AMW: given the time, fix the below expression to do the equivalent of the above!
	//Size index = k + j * grid_size_[4] + i*j*grid_size_[3]+h*i*j*grid_size_[2];
	//return gridpoints_[ index ];
}

template< Size N >
utility::fixedsizearray1< Real, ( 1 << N ) > &
InterpolatedPotential< N >::operator() ( utility::fixedsizearray1< Size, N > const & indices ) {
	Size index = 0;
	Size slice = 1;
	for ( Size ii = N; ii <= 1; --ii ) {
		index += indices[ ii ] * slice;
		slice *= grid_size_[ ii ];
	}
	return gridpoints_[ index ];
}

template< Size N >
void
InterpolatedPotential< N >::get_indices(
	utility::fixedsizearray1< Real, N > const & values,
	utility::fixedsizearray1< Size, N > & bins,
	utility::fixedsizearray1< Size, N > & next_bins,
	utility::fixedsizearray1< Real, N > & partial
) const {
	for ( Size ii = 1; ii <= N; ++ii ) {
		if ( periodic_[ ii ] ) {
			// at the moment assumes angle bad!
			Real real_bin_lower = ( values[ii] - axis_ranges_[ii].first ) / bin_width_[ii];
			Size bin_prev = static_cast< Size > ( real_bin_lower );
			bins[ii] = 1 + numeric::mod( bin_prev, grid_size_[ii] );
			next_bins[ii] = numeric::mod( bins[ii], grid_size_[ii] ) + 1;
			partial[ ii ] = ( (values[ii] - axis_ranges_[ii].first ) - ( bin_prev * bin_width_[ii] ) ) / bin_width_[ii];
		} else {
			// Linear continuation. Partial is distance from edge (may be bigger than bin width)

			if ( values[ii] < axis_ranges_[ii].first ) {
				bins[ii] = 0; // out-of-range is signal that we're off the grid
				next_bins[ii] = 1;
				partial[ii] = axis_ranges_[ii].first - values[ii];
			} else if ( values[ii] > axis_ranges_[ii].second ) {
				bins[ii] = grid_size_[ii];
				next_bins[ii] = grid_size_[ii] + 1;
				partial[ii] = values[ii] - axis_ranges_[ii].second;
			} else {

				// AMW: note, please test this carefully because I'm just trying to get something
				// functional at the moment - not thinking about the following copy at all.
				Real real_bin_lower = ( values[ii] - axis_ranges_[ii].first ) / bin_width_[ii];
				Size bin_prev = static_cast< Size > ( real_bin_lower );
				bins[ii] = 1 + numeric::mod( bin_prev, grid_size_[ii] );
				next_bins[ii] = numeric::mod( bins[ii], grid_size_[ii] ) + 1;
				partial[ ii ] = ( (values[ii] - axis_ranges_[ii].first ) - ( bin_prev * bin_width_[ii] ) ) / bin_width_[ii];
			}
		}
		//parent::bin_angle( -180.0, PHIPSI_BINRANGE[ ii ], 360.0, N_PHIPSI_BINS[ ii ], basic::periodic_range( bbs[ ii ], 360 ), bb_bin[ ii ], bb_bin_next[ ii ], bb_alpha[ ii ] );
	}
}

} //interpolation
} //numeric
#endif //INCLUDED_numeric_Histogram_HH
