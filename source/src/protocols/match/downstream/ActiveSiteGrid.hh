// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/ActiveSiteGrid.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_downstream_ActiveSiteGrid_hh
#define INCLUDED_protocols_match_downstream_ActiveSiteGrid_hh

// Unit headers

// Package headers
#include <protocols/match/BumpGrid.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

namespace protocols {
namespace match {
namespace downstream {

class ActiveSiteGrid : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::ReferenceCount parent;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;

public:
	virtual ~ActiveSiteGrid();

	ActiveSiteGrid();
	ActiveSiteGrid( ActiveSiteGrid const & );

	ActiveSiteGrid const &
	operator = ( ActiveSiteGrid const & rhs );

	/// @brief Set the bounding box for this grid
	void
	set_bounding_box( BoundingBox const & bb );

	void
	set_bin_width( Real width );

	/// @brief Accessor for the bounding box
	BoundingBox const &
	bounding_box() const {
		return bb_;
	}

	Bool3DGrid const &
	grid() const;

	/// @brief Is a point in this grid active?  False for a point outside the bounding box.
	bool
	occupied( Vector const & p ) const;

	/// @brief Reset all the voxels to false
	void
	clear();

	void
	initialize_from_gridlig_file( std::string const & fname );

	/// @brief Set the bounding box to be large enough to hold the volume within the
	/// radius of any atom in the given residue.  This function has the side-effect of
	/// clearing the grid.
	void
	enlargen_to_capture_volume_within_radius_of_residue(
		core::conformation::Residue const & res,
		Real radius
	);

	/// @brief Set the bounding box to be large enough to hold the volume within the
	/// radius of any sidechain atom in the given residue.  This function has the side-effect of
	/// clearing the grid.
	void
	enlargen_to_capture_volume_within_radius_of_sidechain(
		core::conformation::Residue const & res,
		Real radius
	);

	/// @brief Set the bounding box to be large enough to hold the volume within the
	/// radius of any backbone atom in the given residue.  This function has the side-effect of
	/// clearing the grid.
	void
	enlargen_to_capture_volume_within_radius_of_backbone(
		core::conformation::Residue const & res,
		Real radius
	);


	/// @brief Set all the voxels within a certain radius of the residue atoms to true.
	void
	or_within_radius_of_residue(
		core::conformation::Residue const & res,
		Real radius
	);

	/// @brief Set all the voxels within a certain radius of the sidechain atoms to true.
	void
	or_within_radius_of_sidechain(
		core::conformation::Residue const & res,
		Real radius
	);

	/// @brief Set all the voxels within a certain radius of the backbone atoms to true.
	void
	or_within_radius_of_backbone(
		core::conformation::Residue const & res,
		Real radius
	);

	/// @brief Ensures the grid is up-to-date after any calls to enlargen_*.
	void initialize();

private:

	/// @brief Increase the size of the bounding box to contain a given volume
	void
	swell_bounding_box(
		Vector lower,
		Vector upper,
		Real radius
	);

	/// @brief Ensure that the grid has the appropriate bounding box.  This resets the grid
	/// if the object has modifications to its bounding box that have not yet been transfered
	/// to its grid.
	void
	prep_grid();

private:

	Real bin_width_;
	BoundingBox bb_;
	Bool3DGridOP grid_;
	bool reset_grid_bb_;

};

}
}
}

#endif
