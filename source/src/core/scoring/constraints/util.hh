// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/util.hh
/// @brief Utility functions for defining and using constraints.
/// @author James Thompson
/// @author Steven Lewis smlewi@gmail.com (merge_constraints_from_cmdline...)

#ifndef INCLUDED_core_scoring_constraints_util_hh
#define INCLUDED_core_scoring_constraints_util_hh

#include <core/types.hh>

#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>

#include <utility/vector1.hh>
#include <iostream>

#ifdef WIN32
#include <string>
#endif


namespace core {
namespace scoring {
namespace constraints {

//////////////////MATH///////////////////////////////////

/// @brief Returns the log of the weighted value of a Gaussian distribution
/// evaluated with the given mean, sd, and x values. Returns zero if the
/// weight is less than 1e-10.
Real logdgaussian_deriv( Real x, Real mean, Real sd, Real weight );

/// @brief Returns the log of the weighted value of a Gaussian distribution
/// evaluated with the given mean, sd, and x values. Returns zero if the
/// weight is less than 1e-10.
Real logdgaussian( Real x, Real mean, Real sd, Real weight );

// @brief Returns the weighted value of a Gaussian distribution evaluated
// with the given mean, sd, and x values. Returns zero if the weight is less
// than 1e-10.
Real dgaussian( Real x, Real mean, Real sd, Real weight );

/// @brief Returns the weighted derivative of a Gaussian distribution
/// evaluated with the given mean, sd, and x values. Returns zero if the
/// weight is less than 1e-10.
Real gaussian_deriv( Real x, Real mean, Real sd, Real weight );

/// @brief Returns the weighted value of an Exponential distribution evaluated
/// with the given anchor, rate, and x values. Returns zero if the weight is
/// less than 1e-10.
Real dexponential( Real x, Real anchor, Real rate, Real weight );

/// @brief Returns the weighted derivative of an Exponential distribution
/// evaluated with the given anchor, rate, and x values. Returns zero if the
/// weight is less than 1e-10.
Real exponential_deriv( Real x, Real anchor, Real rate, Real weight );

/// @brief Estimates the y-value of the given x-value by interpolating between
/// the given points (x1,y1) and (x2,y2) by using linear interpolation between
/// the two points.
Real linear_interpolate(
	Real const x_val,
	Real const x1,
	Real const x2,
	Real const y1,
	Real const y2
);

//return only those constraints that evaluate with less than threshold on
//the filter_pose
void
cull_violators(
	ConstraintCOPs const & target_list,
	ConstraintCOPs & culled_list,
	core::pose::Pose const & filter_pose,
	core::Real threshold = 1.0
);

///////////////COMMAND LINE/////////////////////////////

////////// Centroid constraints (add and replace)

std::string get_cst_file_option();
//// @brief  add constraints if specified by user.
void add_constraints_from_cmdline_to_pose( core::pose::Pose & pose );
//// @brief  add constraints if specified by user.
void add_constraints_from_cmdline_to_scorefxn(
	core::scoring::ScoreFunction & scorefxn_
);
//// @brief  add constraints if specified by user.
void add_constraints_from_cmdline(
	core::pose::Pose & pose, core::scoring::ScoreFunction & scorefxn_
);

////////// FA constraints (add and replace)

std::string get_cst_fa_file_option();

/// @brief add constraints if specified by user.
void add_fa_constraints_from_cmdline_to_pose( core::pose::Pose & pose );

/// @brief add constraints if specified by user.
void add_fa_constraints_from_cmdline_to_scorefxn(
	core::scoring::ScoreFunction & scorefxn_
);

/// @brief add constraints if specified by user.
void add_fa_constraints_from_cmdline(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn_
);


////////// Centroid constraints (merge mode)

/// @brief merge cmdline constraints to pre-existing constraints
void merge_constraints_from_cmdline_to_pose( core::pose::Pose & pose );

/// @brief merge cmdline constraints to pre-existing constraints - only adds to ZERO weights; previously nonzero constraint weights are unmodified and a warning is issued
void merge_constraints_from_cmdline_to_scorefxn(
	core::scoring::ScoreFunction & scorefxn_
);
/// @brief merge cmdline constraints to pre-existing constraints
void merge_constraints_from_cmdline(
	core::pose::Pose & pose, core::scoring::ScoreFunction & scorefxn_
);

////////// FA constraints (merge mode)
/// @brief merge cmdline constraints to pre-existing constraints
void merge_fa_constraints_from_cmdline_to_pose( core::pose::Pose & pose );

/// @brief merge cmdline constraints to pre-existing constraints - only adds to ZERO weights; previously nonzero constraint weights are unmodified and a warning is issued
void merge_fa_constraints_from_cmdline_to_scorefxn(
	core::scoring::ScoreFunction & scorefxn_
);

/// @brief merge cmdline constraints to pre-existing constraints
void merge_fa_constraints_from_cmdline(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn_
);


/////////////////Coordinate Constraints

/// @brief handy function for tethering pose to starting coordinates.
void
add_coordinate_constraints( core::pose::Pose & pose, core::Real const coord_sdev = 10.0, bool include_sc = true);

/// @brief Add coordinate constraints for starting coordinates, start:end residues, inclusive.
void
add_coordinate_constraints( core::pose::Pose & pose, core::Size const start_res, core::Size const end_res, core::Real const coord_sdev = 10.0, bool include_sc = true);


////////// Constraint removal

/// @brief Remove all constraints of a given type from a pose.
void
remove_constraints_of_type(core::pose::Pose & pose, std::string const & type);

/// @brief Remove all constraints of a given type from a pose that involve start_res to end_res.  Useful for coordinate/dihedral constraints
void
remove_constraints_of_type(core::pose::Pose & pose, std::string const & type, core::Size const start_res, core::Size const end_res);

void
remove_nonbb_constraints( pose::Pose & pose) ;


/// @brief call this on your constraints if you have MultiConstraints before running Abinitio -- already done by broker-type application
void choose_effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp, ConstraintCOPs& in );

/// @brief combine constraints randomly into Ambiguous constraints...
void combine_constraints(
	ConstraintCOPs& in,
	core::Size combine_ratio,
	utility::vector1< bool > exclude_res,
	core::kinematics::ShortestPathInFoldTree const& sp
);

/// @brief have at most one constraint per residue pair...
void skip_redundant_constraints( ConstraintCOPs& in, core::Size total_residue, core::Size influence_width = 1 );
void drop_constraints( ConstraintCOPs& in, core::Real drop_rate );

/// @brief example of how to go through a pose constraint set and print out stuff.
void print_atom_pair_constraints( pose::Pose const & pose, std::ostream & out = std::cout );


/// @brief map constraints to new atom numbers after, e.g. variants change. requires pose to have same number of residues.
void
map_constraints_from_original_pose( core::pose::Pose const & original_pose, core::pose::Pose & pose );


} // namespace constraints
} // namespace scoring
} // namespace core

#endif
