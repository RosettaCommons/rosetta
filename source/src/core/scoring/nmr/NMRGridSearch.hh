// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRGridSearch.hh
/// @brief   class that performs a grid search to find the metal ion position that is in best
///          agreement with PCS and PRE data
/// @details last Modified: 06/16/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRGridSearch_HH
#define INCLUDED_core_scoring_nmr_NMRGridSearch_HH

// Unit headers
#include <core/scoring/nmr/NMRGridSearch.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

namespace core {
namespace scoring {
namespace nmr {

class NMRGridSearch {

public: // Methods

	/// @brief constructor which takes as arguments two AtomIDs
	///        that define the location of the grid search.
	///        the area in which the grid search is performed is defined by
	///        default values for the distance_to_atom1 (= 5 Ang), the stepsize (= 3 Ang),
	///        the grid_min_radius (= 2.5 Ang) and the grid_max_radius (= 25 Ang).
	///        to define a different area use the public set methods.
	NMRGridSearch(
		id::AtomID const & grid_atom1,
		id::AtomID const & grid_atom2,
		pose::Pose const & pose
	);

	/// @brief constructor with a full list of arguments that define the
	///        location and area of the grid search.
	NMRGridSearch(
		id::AtomID const & grid_atom1,
		id::AtomID const & grid_atom2,
		pose::Pose const & pose,
		Real const distance_center_to_atom1,
		Real const stepsize,
		Real const grid_min_radius,
		Real const grid_max_radius
	);

	/// @brief copy constructor
	NMRGridSearch(NMRGridSearch const & other);

	/// @brief assignment operator
	NMRGridSearch&
	operator=(NMRGridSearch const & rhs);

	/// @ destructor
	~NMRGridSearch();

	/// @brief checks if grid point is within the predefined sphere and updates metal coordinates
	///        returns false when last grid point is reached
	bool
	valid_next_grid_point(Vector & metal_coords);

	// Getters
	id::AtomID const & get_grid_atom1() const { return grid_atom1_; }
	id::AtomID const & get_grid_atom2() const { return grid_atom2_; }
	Real get_stepsize() const { return stepsize_; }
	Real get_grid_min_radius() const { return grid_min_radius_; }
	Real get_grid_max_radius() const { return grid_max_radius_; }
	Real get_distance_center_to_atom1() const { return distance_center_to_atom1_; }
	PointPosition get_grid_search_center() const { return center_; }
	PointPosition get_best_grid_point() const { return best_grid_point_; }

	// Setters
	void set_stepsize(Real step) { stepsize_ = step; }
	void set_grid_min_radius(Real radius) { grid_min_radius_ = radius; }
	void set_grid_max_radius(Real radius) { grid_max_radius_ = radius; }
	void set_distance_center_to_atom1(Real distance) { distance_center_to_atom1_ = distance; }
	void set_grid_atom1(id::AtomID const & atom) { grid_atom1_ = atom; }
	void set_grid_atom2(id::AtomID const & atom) { grid_atom2_ = atom; }
	void set_best_grid_point(PointPosition const & point) { best_grid_point_ = point; }

	/// @brief Calculate and set grid search center from coordinates
	///        of atoms 1 and 2 found in the input pose and from
	///        the geometric definition of the grid search
	void set_grid_search_center(pose::Pose const & pose);

	/// @brief Set grid search center directly to specified input point
	void set_grid_search_center(PointPosition const & point);

private: // Methods

	/// @brief default constructor
	NMRGridSearch();

	void reset_to_start();

	/// @brief advances metal coordinates by one stepsize in either x, y or z-direction
	/// returns false when the last grid point is reached
	bool switch_to_next_grid_point(Vector & metal_coords);

private: // Data

	id::AtomID grid_atom1_;
	id::AtomID grid_atom2_;
	Real distance_center_to_atom1_;
	PointPosition center_;
	PointPosition best_grid_point_; // grid point with lowest score
	PointPosition current_;
	Real stepsize_;
	Real delta_;
	Real grid_min_radius_;
	Real grid_max_radius_;
	bool center_is_set_;
	bool at_start_position_;

};

} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRGridSearch_HH
