// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRGridSearch.cc
/// @brief   Implementation of class NMRGridSearch
/// @details last Modified: 06/21/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/NMRGridSearch.hh>

// Package headers
#include <core/io/nmr/AtomSelection.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <cmath>

namespace core {
namespace scoring {
namespace nmr {

static basic::Tracer TR( "core.scoring.nmr.NMRGridSearch" );

/// @brief constructor which takes as arguments two AtomIDs
///        that define the location of the grid search.
///        the area in which the grid search is performed is defined by
///        default values for the distance_to_atom1 (= 5 Ang), the stepsize (= 3 Ang),
///        the grid_min_radius (= 2.5 Ang) and the grid_max_radius (= 25 Ang).
///        to define a different area use the public set methods.
NMRGridSearch::NMRGridSearch(
	id::AtomID const & grid_atom1,
	id::AtomID const & grid_atom2,
	pose::Pose const & pose
) :
	grid_atom1_(grid_atom1),
	grid_atom2_(grid_atom2),
	distance_center_to_atom1_(5.0),
	center_(),
	best_grid_point_(),
	stepsize_(3.0),
	delta_(1.5),
	grid_min_radius_(2.5),
	grid_max_radius_(25.0),
	center_is_set_(false),
	at_start_position_(false)
{
	// We want to be able to set the grid search centered on one atom
	//runtime_assert_msg( grid_atom1_ != grid_atom2_ , "ERROR: Grid atom 1 and 2 must be different." );
	set_grid_search_center(pose);
}

/// @brief constructor with a full list of arguments that define the
///        location and area of the grid search.
NMRGridSearch::NMRGridSearch(
	id::AtomID const & grid_atom1,
	id::AtomID const & grid_atom2,
	pose::Pose const & pose,
	Real const distance_center_to_atom1,
	Real const stepsize,
	Real const grid_min_radius,
	Real const grid_max_radius
) :
	grid_atom1_(grid_atom1),
	grid_atom2_(grid_atom2),
	distance_center_to_atom1_(distance_center_to_atom1),
	center_(),
	best_grid_point_(),
	stepsize_(stepsize),
	delta_(stepsize/2.0),
	grid_min_radius_(grid_min_radius),
	grid_max_radius_(grid_max_radius),
	center_is_set_(false),
	at_start_position_(false)
{
	runtime_assert_msg( stepsize > 0.0, "ERROR: Step size in NMRGridSearch must be positive." );
	runtime_assert_msg( grid_max_radius > grid_min_radius, "ERROR: Outer sphere radius must be bigger than inner sphere radius." );
	//runtime_assert_msg( grid_atom1_ != grid_atom2_, "ERROR: Grid atom 1 and 2 must be different." );
	set_grid_search_center(pose);
}

/// @brief copy constructor
NMRGridSearch::NMRGridSearch(NMRGridSearch const & other) :
	grid_atom1_(other.grid_atom1_),
	grid_atom2_(other.grid_atom2_),
	distance_center_to_atom1_(other.distance_center_to_atom1_),
	center_(other.center_),
	best_grid_point_(other.best_grid_point_),
	current_(other.current_),
	stepsize_(other.stepsize_),
	delta_(other.delta_),
	grid_min_radius_(other.grid_min_radius_),
	grid_max_radius_(other.grid_max_radius_),
	center_is_set_(other.center_is_set_),
	at_start_position_(other.at_start_position_)
{}

/// @brief assignment operator
NMRGridSearch&
NMRGridSearch::operator=(NMRGridSearch const & rhs) {
	if ( this != &rhs ) {
		grid_atom1_ = rhs.grid_atom1_;
		grid_atom2_ = rhs.grid_atom2_;
		distance_center_to_atom1_ = rhs.distance_center_to_atom1_;
		center_ = rhs.center_;
		best_grid_point_ = rhs.best_grid_point_;
		current_ = rhs.current_;
		stepsize_ = rhs.stepsize_;
		delta_ = rhs.delta_;
		grid_min_radius_ = rhs.grid_min_radius_;
		grid_max_radius_ = rhs.grid_max_radius_;
		center_is_set_ = rhs.center_is_set_;
		at_start_position_ = rhs.at_start_position_;
	}
	return *this;
}

/// @ destructor
NMRGridSearch::~NMRGridSearch() { }

void
NMRGridSearch::reset_to_start() {
	current_ = center_;
	at_start_position_ = true;
}

/// @brief checks if grid point is within the predefined sphere and updates metal coordinates
///        returns false when last grid point is reached
bool
NMRGridSearch::valid_next_grid_point(Vector & metal_coords) {
	if ( !center_is_set_ ) {
		utility_exit_with_message("ERROR: NMRGridSearch center is not set. Call \"set_grid_search_center()\" before starting grid search.");
	}
	while ( switch_to_next_grid_point(metal_coords) == true ) {
		Real r = metal_coords.distance(center_);
		if ( ( r <= grid_max_radius_ + delta_ ) && ( r >= grid_min_radius_ - delta_ ) ) {
			return true;
		}
	}
	return false;
}

/// @brief advances metal coordinates by one stepsize in either x, y or z-direction
///        returns false when the last grid point is reached
bool
NMRGridSearch::switch_to_next_grid_point(Vector & metal_coords) {
	metal_coords = current_;

	if ( at_start_position_ ) {
		current_ = center_ - grid_max_radius_;
		at_start_position_ = false;

		// Any move would place the point outside of the box
		if ( stepsize_ > (2.0*grid_max_radius_ + delta_) ) {
			reset_to_start();
			return false;
		}
		return true;
	}

	// Advance point by one step in x-direction.
	current_.x() += stepsize_;

	// Max allowed x-coordinate is reached.
	// Move back to min allowed grid point in x-direction
	// and advance grid point by one step in y-direction.
	if ( (current_.x() - center_.x() ) > grid_max_radius_ + delta_ ) {
		current_.x( -grid_max_radius_ + center_.x() );
		current_.y() += stepsize_;
	}
	// Max allowed y-coordinate is reached.
	// Move back to min allowed grid point in y-direction
	// and advance grid point by one step in z-direction.
	if ( (current_.y() - center_.y() ) > grid_max_radius_ + delta_ ) {
		current_.y( -grid_max_radius_ + center_.y() );
		current_.z() += stepsize_;
	}
	// Grid search has reached end and point cannot
	// be advanced any further. So it is not a valid next
	// grid point. Reset grid search to start.
	if ( (current_.z() - center_.z() ) > grid_max_radius_ + delta_ ) {
		current_.z( -grid_max_radius_ + center_.z() );
		reset_to_start();
		return false;
	}
	return true;
}
/// @brief Calculate and set grid search center from coordinates
///        of atoms 1 and 2 found in the input pose and from
///        the geometric definition of the grid search
void
NMRGridSearch::set_grid_search_center(pose::Pose const & pose) {
	if ( !grid_atom1_.valid() || grid_atom1_.rsd() > pose.total_residue() ) {
		utility_exit_with_message( "ERROR in updating grid search center. First grid search atom not found. Check NMR input file." );
	}
	if ( !grid_atom2_.valid() || grid_atom2_.rsd() > pose.total_residue() ) {
		utility_exit_with_message( "ERROR in updating grid search center. Second grid search atom not found. Check NMR input file." );
	}

	PointPosition atom1_coords(pose.residue(grid_atom1_.rsd()).xyz(grid_atom1_.atomno()));
	PointPosition atom2_coords(pose.residue(grid_atom2_.rsd()).xyz(grid_atom2_.atomno()));
	if ( grid_atom1_ == grid_atom2_ ) {
		center_ = atom1_coords;
	} else {
		center_ = atom1_coords + (distance_center_to_atom1_/atom1_coords.distance(atom2_coords)) * (atom2_coords - atom1_coords);
	}
	center_is_set_ = true;
	reset_to_start();
}

/// @brief Set grid search center directly to specified input point
void
NMRGridSearch::set_grid_search_center(PointPosition const & point)
{
	center_ = point;
	center_is_set_ = true;
	reset_to_start();
}

} // namespace nmr
} // namespace scoring
} // namespace core
