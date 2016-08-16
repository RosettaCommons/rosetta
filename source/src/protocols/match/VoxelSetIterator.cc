// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/VoxelSetIterator.cc
/// @brief  Implementation of the VoxelSetIterator class
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/VoxelSetIterator.hh>

// Utility headers
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>

// C++ headers

namespace protocols {
namespace match {

using namespace utility;

/// @details Initalize the iterator with both the bounding volume and discritezation data
/// from the OccupiedSpaceHasher, and also with the coordinates of the 6-d point that this
/// instance will be responsible for.  One point per VoxelSetIterator.
VoxelSetIterator::VoxelSetIterator(
	BoundingBox const & bb,
	Size3 const & n_xyz_bins,
	Size3 const & n_euler_bins,
	Real3 const & xyz_bin_widths,
	Real3 const & euler_bin_widths,
	Real3 const & xyz_bin_halfwidths,
	Real3 const & euler_bin_halfwidths,
	Real6 const & point
) :
	bb_( bb ),
	n_xyz_bins_( n_xyz_bins ),
	n_euler_bins_( n_euler_bins ),
	xyz_bin_widths_( xyz_bin_widths ),
	euler_bin_widths_( euler_bin_widths ),
	xyz_bin_halfwidths_( xyz_bin_halfwidths ),
	euler_bin_halfwidths_( euler_bin_halfwidths ),
	point_( point ),
	//wrap_theta_( false ),
	theta_near_0_( false ),
	theta_near_180_( false ),
	wrapped_phipsi_bins_( /* 0 */ ),
	wrapped_phipsi_halfbins_( /* 0 */ ),
	curr_bin_( /* 0 */ ),
	curr_pos_( 0 )
{
	assert( point_[ 4 ] >= 0 && point_[ 4 ] <  360 );
	assert( point_[ 5 ] >= 0 && point_[ 5 ] <  360 );
	assert( point_[ 6 ] >= 0 && point_[ 6 ] <= 180 );

	/// OK.  So we're going to enumerate all 64 voxels that this point will hash into.
	/// This is made slightly tricky by three facts:
	/// 1. We might walk outside of the bounding box in the xyz dimensions
	/// 2. For euler angles "phi" and "psi", the angles are periodic in the range between [ 0 .. 360 ]
	/// 3. For euler angle "theta" near 180 or 0, we "flip" the "phi" and "psi" angles if phi is greater than 180.

	for ( Size ii = 1; ii <= 3; ++ii ) {
		basebin_[ ii ]     = static_cast< Size > (( point_[ ii ] - bb_.lower()( ii )) / xyz_bin_widths_[ ii ] );
		basehalfbin_[ ii ] = static_cast< Size > (( point_[ ii ] - bb_.lower()( ii )) / xyz_bin_halfwidths_[ ii ] ) - 2 * basebin_[ ii ];
		if ( basebin_[ ii ] >= n_xyz_bins_[ ii ] ) return;
	}
	for ( Size ii = 4, ii_minus3 = 1; ii <= 6; ++ii, ++ii_minus3 ) {
		assert( point_[ ii ] >= 0 && point_[ ii ] <= 360 );
		assert( ii != 6 || point_[ ii ] <= 180 ); /// theta must range between 0 and 180.
		basebin_[ ii ]     = static_cast< Size > ( point_[ ii ] / euler_bin_widths_[ ii_minus3 ] );
		basehalfbin_[ ii ] = static_cast< Size > ( point_[ ii ] / euler_bin_halfwidths_[ ii_minus3 ] ) - 2 * basebin_[ ii ];

		// In floating point, sometimes (360 - epsilon) / X = 360 / X for
		// values of X that evenly divide 360, and epsilon non-zero. )
		if ( ii != 6 && basebin_[ ii ] == n_euler_bins_[ ii_minus3 ] ) {
			basebin_[ ii ] = n_euler_bins_[ ii_minus3 ] - 1;
		}

	}

	fixedsizearray1< Size, 6 > twos( 2 );
	iter64_.set_dimension_sizes( twos );

	theta_near_0_   = point_[ 6 ] < euler_bin_halfwidths[ 3 ];
	theta_near_180_ = point_[ 6 ] > 180 - euler_bin_halfwidths[ 3 ];

	if ( theta_near_0_ || theta_near_180_ ) {
		Real wrapped_phi, wrapped_psi;
		if ( point_[ 4 ] > 180 ) {
			/// wrap both phi and psi to their negative-rotation values
			/// when the original phi value is greater than 180; the new
			/// phi value will be less than 180.  The new psi may end up negative, so check if psi - 180 is < 0.
			wrapped_phi = point_[ 4 ] - 180;
			wrapped_psi = point_[ 5 ] > 180 ? point_[ 5 ] - 180 : 180 + point_[ 5 ];
		} else {
			/// keep phi and psi fixed since phi is < 180 -- let the neighbors > 180 come to us.
			wrapped_phi = point_[ 4 ];
			wrapped_psi = point_[ 5 ];
		}
		wrapped_phipsi_bins_[ 1 ]     = static_cast< Size > (( wrapped_phi ) / euler_bin_widths_[ 1 ] );
		wrapped_phipsi_bins_[ 2 ]     = static_cast< Size > (( wrapped_psi ) / euler_bin_widths_[ 2 ] );
		wrapped_phipsi_halfbins_[ 1 ] = static_cast< Size > (( wrapped_phi ) / euler_bin_halfwidths_[ 1 ] ) - 2 * wrapped_phipsi_bins_[ 1 ];
		wrapped_phipsi_halfbins_[ 2 ] = static_cast< Size > (( wrapped_psi ) / euler_bin_halfwidths_[ 2 ] ) - 2 * wrapped_phipsi_bins_[ 2 ];
	}
	calc_bin_and_pos();
}

void VoxelSetIterator::operator ++ ()
{
	if ( iter64_.at_end() ) return;
	Size start = 7 - ++iter64_; /// start at the x coordinate if all 6 dimensions rolled over;
	while ( true ) {
		if ( iter64_.at_end() ) break;

		bool outside_bounding_box( false );
		for ( Size ii = start; ii <= 3; ++ii ) {
			if ( iter64_[ ii ] == 2 && basebin_[ ii ] == n_xyz_bins_[ ii ] - 1 && basehalfbin_[ ii ] == 1 ) {
				/// walked out of bounds.
				start = 7 - iter64_.continue_at_dimension( ii );
				outside_bounding_box = true;
				break;
			}
		}
		if ( ! outside_bounding_box ) break;
	}

	if ( ! iter64_.at_end() ) calc_bin_and_pos();
}

bool VoxelSetIterator::at_end() const
{
	return iter64_.at_end();
}

void VoxelSetIterator::get_bin_and_pos( Size6 & bin, Size & pos ) const
{
	assert( ! at_end() );
	bin = curr_bin_;
	pos = curr_pos_;
}

void VoxelSetIterator::calc_bin_and_pos()
{
	assert( ! at_end() );
	curr_pos_ = 1; /// curr pos ranges from 1 and 64.  The math below will add between 0 and 63 to curr_pos_.

	for ( Size ii = 1; ii <= 3; ++ii ) {
		curr_bin_[ ii ] = basebin_[ ii ];
		if ( iter64_[ ii ] == 2 && basehalfbin_[ ii ] == 1 ) {
			curr_bin_[ ii ] += 1;
		} else if ( basehalfbin_[ ii ] == 1 || iter64_[ ii ] == 2 ) {
			curr_pos_ += ii == 1 ? 32 : ( ii == 2 ? 16 : 8 );
		}
	}

	/// Now handle the Euler bins.
	/// Two cases:
	/// 1. Near 0 or Near 360
	/// 2. Anywhere else

	if ( ( theta_near_0_ && iter64_[ 6 ] == 1 ) || ( theta_near_180_ && iter64_[ 6 ] == 2 ) ) {

		// phi and psi
		for ( Size ii = 4, ii_minus3 = 1; ii <= 5; ++ii, ++ii_minus3 ) {
			curr_bin_[ ii ] = wrapped_phipsi_bins_[ ii_minus3 ];
			if ( iter64_[ ii ] == 2 && wrapped_phipsi_halfbins_[ ii_minus3 ] == 1 ) {
				// increment, and wrap at 360 if necessary
				curr_bin_[ ii ] = curr_bin_[ ii ] + 1 == n_euler_bins_[ ii_minus3 ] ? 0 : curr_bin_[ ii ] + 1;
			} else if ( iter64_[ ii ] == 2 || wrapped_phipsi_halfbins_[ ii_minus3 ] == 1 ) {
				curr_pos_ += ii == 4 ? 4 : 2;
			}
		}

		/// theta
		if ( theta_near_0_ ) {
			curr_bin_[ 6 ] = 0;
		} else {
			curr_bin_[ 6 ] = n_euler_bins_[ 3 ] - 1;
			curr_pos_ += 1;
		}
	} else {
		for ( Size ii = 4, ii_minus3 = 1; ii <= 6; ++ii, ++ii_minus3 ) {
			curr_bin_[ ii ] = basebin_[ ii ];
			if ( iter64_[ ii ] == 2 && basehalfbin_[ ii ] == 1 ) {
				// Wrap at > 360 -- only applies to phi and psi, but we've gotten to this point knowing
				// that the above condition does not apply to theta, so no need for an ii==4 || ii==5 check.
				assert( ii != 6 || curr_bin_[ ii ] + 1 != n_euler_bins_[ ii_minus3 ]  );
				curr_bin_[ ii ] = ( curr_bin_[ ii ] + 1 == n_euler_bins_[ ii_minus3 ] ? 0 : curr_bin_[ ii ] + 1 );
			} else if ( basehalfbin_[ ii ] == 1 || iter64_[ ii ] == 2 ) {
				curr_pos_ += ii == 4 ? 4 : ( ii == 5 ? 2 : 1 );
			}
		}
	}
}

}
}
