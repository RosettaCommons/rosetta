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

#ifndef INCLUDED_numeric_InterpolatedPotential_HH
#define INCLUDED_numeric_InterpolatedPotential_HH

// Unit Headers
#include <numeric/interpolation/InterpolatedPotential.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

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
class InterpolatedPotential {

public:

	InterpolatedPotential():
		gridpoints_( nullptr )
	{
		for ( auto & p : periodic_ ) { p = false; }
	}

	~InterpolatedPotential() { delete [] gridpoints_; }

public:
	void dimension( utility::fixedsizearray1< Size, N > const & dims ) {
		bool ret = true;
		for ( Size ii = 1; ii <= N; ++ii ) {
			if ( dims[ ii ] != grid_size_[ ii ] ) ret = false;
		}
		if ( ret ) return; // noop

		Size product = 1;
		for ( Size const dim : dims ) product *= dim;
		if ( gridpoints_ ) delete [] gridpoints_;
		gridpoints_ = new utility::fixedsizearray1< Real, ( 1 << N ) >[ product ];
		grid_size_ = dims;
	}

public:

	utility::fixedsizearray1< Real, ( 1 << N ) > &
	operator() ( Size i, Size j, Size k, Size l );

	utility::fixedsizearray1< Real, ( 1 << N ) > &
	operator() ( utility::fixedsizearray1< Size, N > const & indices );

	utility::fixedsizearray1< Real, ( 1 << N ) > const &
	operator() ( utility::fixedsizearray1< Size, N > const & indices ) const;

public:

	utility::fixedsizearray1< Size, N > const &
	grid_size() const {
		return grid_size_;
	}

	std::pair< Real, Real > const &
	axis_range( Size const n ) const {
		return axis_ranges_[ n ];
	}

	utility::fixedsizearray1< Real, N > const &
	bin_width() const {
		return bin_width_;
	}

	void set_bin_width( utility::fixedsizearray1< Real, N > const & values ) {
		bin_width_ = values;
	}

	void
	get_indices(
		utility::fixedsizearray1< Real, N > const & values,
		utility::fixedsizearray1< Size, N > & bins,
		utility::fixedsizearray1< Size, N > & next_bins,
		utility::fixedsizearray1< Real, N > & partial
	) const;

	void set_periodic( bool const setting ) {
		for ( bool & dim : periodic_ ) dim = setting;
	}

	void set_periodic( utility::fixedsizearray1< bool, N > const & periodic ) { periodic_ = periodic; }

private:
	utility::fixedsizearray1< Size, N > grid_size_;
	utility::fixedsizearray1< std::pair< Real, Real >, N > axis_ranges_ = std::make_pair( 0.0, 0.0 );
	//std::array< utility::fixedsizearray1< Real, ( 1 << N ) > > gridpoints_;
	// Implementing as a raw pointer. We can't know how many gridpoints there are
	// yet. Look into this later, but you'd need variadic templates maybe?
	utility::fixedsizearray1< Real, ( 1 << N ) > * gridpoints_ = nullptr;

	utility::fixedsizearray1< Real, N > bin_width_;
	utility::fixedsizearray1< bool, N > periodic_;

};

template< Size N >
void
polycubic_interpolation(
	numeric::interpolation::InterpolatedPotential< N > const & interpolated_potential,
	utility::fixedsizearray1< Real, N > const & values,
	Real & score,
	utility::fixedsizearray1< Real, N > & dscoreddof ) {

	// What indices are we between?
	utility::fixedsizearray1< Size, N > bins;
	utility::fixedsizearray1< Size, N > next_bins;
	utility::fixedsizearray1< Real, N > partial;
	interpolated_potential.get_indices( values, bins, next_bins, partial );

	// Catch linear continuation + aperiodic case.
	bool do_linear = false;
	for ( Size ii = 1; ii <= N; ++ii ) {
		if ( bins[ii] == 0 ) {
			do_linear = true;
			// Reset to being just on the border!
			bins[ii] = 1;
			next_bins[ii] = 2;
			partial[ii] = 0.0;
		}
		if ( next_bins[ii] == interpolated_potential.grid_size()[ii]+1 ) {
			do_linear = true;
			// Reset to being just on the border!
			bins[ii] = interpolated_potential.grid_size()[ii] - 1;
			next_bins[ii] = interpolated_potential.grid_size()[ii];
			partial[ii] = 1.0;
		}
	}

	// construct n_derivs
	utility::fixedsizearray1< utility::fixedsizearray1< Real, ( 1 << N ) >, ( 1 << N ) > n_derivs;
	for ( Size ii = 1; ii <= ( 1 << N ); ++ ii ) {

		utility::fixedsizearray1< Size, N > indices;
		for ( Size jj = 1; jj <= N; ++ii ) {
			// does indices take from this or next
			indices[ jj ] = ( ( ii - 1 ) & ( 1 << ( N - jj ) ) ) ? next_bins[ jj ] : bins[ jj ];
		}
		n_derivs[ ii ] = interpolated_potential( indices );
	}

	// dispatch to polycubic interpolation function
	numeric::interpolation::polycubic_interpolation( n_derivs, partial, interpolated_potential.bin_width(), score, dscoreddof );

	if ( do_linear ) {
		// find which dimensions are out of bounds, and linearly interpolate over each one.
		for ( Size ii = 1; ii <= N; ++ii ) {
			if ( values[ii] < interpolated_potential.axis_range( ii ).first ) {
				score += ( values[ii] - interpolated_potential.axis_range( ii ).first ) * dscoreddof[ ii ];
			} else if ( values[ii] > interpolated_potential.axis_range( ii ).second ) {
				score += ( values[ii] - interpolated_potential.axis_range( ii ).second ) * dscoreddof[ ii ];
			}
		}
	}
}


} //interpolation
} //numeric
#endif //INCLUDED_numeric_InterpolatedPotential_HH
