// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/relax/cst_util.hh
/// @brief small bundle of utilities for incorporating constraints into relax
/// @author James Thompson

#ifndef INCLUDED_protocols_relax_cst_util_hh
#define INCLUDED_protocols_relax_cst_util_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace relax {

void coordinate_constrain_selection(
	core::pose::Pose & pose,
	core::sequence::SequenceAlignment aln,
	core::Real coord_sdev
);

/// @brief Generate a set of coordinate constraints to backbone atoms using the
/// given standard deviations, with one sd per-reside. If no constraint should
/// be applied to a given residue, give a -1 for the value of the sdev.
core::scoring::constraints::ConstraintSetOP
generate_bb_coordinate_constraints(
	core::pose::Pose & pose,
	utility::vector1< core::Real > const & coord_sdevs
);

utility::vector1< core::Real >
get_per_residue_scores(
	core::pose::Pose & pose,
	core::scoring::ScoreType scoretype
);

void add_virtual_residue_to_cterm(
	core::pose::Pose & pose
);

void delete_virtual_residues(
	core::pose::Pose & pose
);

void derive_sc_sc_restraints(
	core::pose::Pose & pose,
	core::Real const upper_dist_cutoff
);

} // relax
} // protocols

#endif
