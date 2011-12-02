// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loops_util_hh
#define INCLUDED_protocols_loops_util_hh

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Package headers
#include <protocols/loops/Loops.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace loops {

//@brief adds coord constraints for the atoms that are not in the loops structure
void fix_with_coord_cst( loops::Loops const& rigid, core::pose::Pose& pose, bool bCstAllAtom, utility::vector1< core::Real >& );

///@brief get frags that are fully within the Loop --- shorten(=true/false) frags that are close to the end of loops.
extern void select_loop_frags(
	loops::Loops const& loops,
	core::fragment::FragSet& source,
	core::fragment::FragSet& loop_frags,
	core::Size min_size = 1 /* set to 0 if you don't want to shorten at all */
);

void set_extended_torsions_and_idealize_loops( core::pose::Pose& pose, loops::Loops loops );

/// @brief Identical to set_extended_torsions_and_idealize_loops() without the irrational
/// behavior surrounding empty loops.
void safe_set_extended_torsions_and_idealize_loops(const protocols::loops::Loops& loops,
                                                   core::pose::Pose* pose);

void addScoresForLoopParts(
	core::pose::Pose& pose,
	loops::Loops loops,
	const core::scoring::ScoreFunction &scorefxn,
	core::pose::Pose& native_pose, core::Size nloops
);

loops::Loops compute_ss_regions(
	core::Real max_loop_frac,
	core::Size min_length,
	core::fragment::SecondaryStructure const & ss
);

core::scoring::ScoreFunctionOP get_fa_scorefxn();
core::scoring::ScoreFunctionOP get_cen_scorefxn();

void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_target_pose,  protocols::loops::Loops &exclude_regions );

/// loop_str has the format: start:end:cut,start:end:cut and can use rosetta or pdb numbering. The return value is a Loops object encoding that loop
Loops
loops_from_string( std::string const loop_str, core::pose::Pose const & pose );

// this function will return a bunch of "loops" that refer to residues that are considered part of the core:
// not scored are loops with 4 or more residues, short helices (<=5) that terminate a loop are not scored, too
void define_scorable_core_from_secondary_structure( core::fragment::SecondaryStructure const&, protocols::loops::Loops& score_core );


} //loops
} //protocols

#endif
