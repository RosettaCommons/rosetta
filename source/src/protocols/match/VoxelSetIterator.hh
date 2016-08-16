// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/VoxelSetIterator.hh
/// @brief  Declaration for iterator to traverse the 64 bins each hit covers in the OccupiedSpaceHash
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_VoxelSetIterator_hh
#define INCLUDED_protocols_match_VoxelSetIterator_hh

// Unit headers
#include <protocols/match/VoxelSetIterator.fwd.hh>

// Package headers
//#include <protocols/match/SixDHasher.fwd.hh>
#include <protocols/match/Hit.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>

// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/FixedSizeLexicographicalIterator.hh>

namespace protocols {
namespace match {

/// @brief Helper class for the OccupiedSpaceHasher which manages the logic for how
/// to iterate across the 64 voxels that each 6-D point covers.
///
/// @details This class ensures that the bounding box for the hash is not
/// walked outside of, that the phi and psi are wrapped at 360, and that
/// when theta is near a gimbal-lock angle of 180 or 0, that phi and psi are
/// appropriately wrapped to the negative rotation.  This class may be rapidly
/// allocated and deallocated on the stack -- no calls to new are made anywhere.
class VoxelSetIterator
{
public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef numeric::geometry::hashing::Real3        Real3;
	typedef numeric::geometry::hashing::Size6        Size6;
	typedef numeric::geometry::hashing::Bin6D        Bin6D;

public:
	VoxelSetIterator(
		BoundingBox const & bb,
		Size3 const & n_xyz_bins,
		Size3 const & n_euler_bins,
		Real3 const & xyz_bin_widths,
		Real3 const & euler_bin_widths,
		Real3 const & xyz_bin_halfwidths,
		Real3 const & euler_bin_halfwidths,
		Real6 const & point
	);

	void operator ++ ();
	bool at_end() const;
	void get_bin_and_pos( Size6 & bin, Size & pos ) const;

private:
	void calc_bin_and_pos();

private:
	BoundingBox bb_;
	Size3 n_xyz_bins_;
	Size3 n_euler_bins_;
	Real3 xyz_bin_widths_;
	Real3 euler_bin_widths_;
	Real3 xyz_bin_halfwidths_;
	Real3 euler_bin_halfwidths_;
	Real6 point_;

	Bin6D basebin_;
	Bin6D basehalfbin_;

	//bool wrap_theta_;
	bool theta_near_0_;
	bool theta_near_180_;

	utility::fixedsizearray1< Size, 2 > wrapped_phipsi_bins_;
	utility::fixedsizearray1< Size, 2 > wrapped_phipsi_halfbins_;

	utility::FixedSizeLexicographicalIterator< 6 > iter64_;
	Size6 curr_bin_; // The voxel for the
	Size  curr_pos_; // which of the 64 different half-width voxels in the regular width voxel
};


}
}

#endif
