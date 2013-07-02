// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/util.hh
/// @brief utility functions for defining constraints. Maybe better placed in src/numeric?
/// @author James Thompson

#ifndef INCLUDED_core_scoring_constraints_util_hh
#define INCLUDED_core_scoring_constraints_util_hh

#include <core/types.hh>

#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>

#include <utility/vector1.hh>


#ifdef WIN32
#include <string>
#endif


namespace core {
namespace scoring {
namespace constraints {

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

////////// Centroid constraints

std::string get_cst_file_option();
//// @brief 	add constraints if specified by user.
void add_constraints_from_cmdline_to_pose( core::pose::Pose & pose );
//// @brief 	add constraints if specified by user.
void add_constraints_from_cmdline_to_scorefxn(
	core::scoring::ScoreFunction & scorefxn_
);
//// @brief 	add constraints if specified by user.
void add_constraints_from_cmdline(
	core::pose::Pose & pose, core::scoring::ScoreFunction & scorefxn_
);

////////// FA constraints

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

/// @brief	handy function for tethering pose to starting coordinates.
void
add_coordinate_constraints( core::pose::Pose & pose, core::Real const coord_sdev = 10.0 );

/// @brief call this on your constraints if you have MultiConstraints before running Abinitio -- already done by broker-type application
void choose_effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp, ConstraintCOPs& in );

///@brief combine constraints randomly into Ambiguous constraints...
void combine_constraints(
  ConstraintCOPs& in,
	core::Size combine_ratio,
	utility::vector1< bool > exclude_res,
	core::kinematics::ShortestPathInFoldTree const& sp
);

///@brief have at most one constraint per residue pair...
void skip_redundant_constraints( ConstraintCOPs& in, core::Size total_residue, core::Size influence_width = 1 );
void drop_constraints( ConstraintCOPs& in, core::Real drop_rate );
void remove_nonbb_constraints( pose::Pose & pose) ;

} // namespace constraints
} // namespace scoring
} // namespace core

#endif
