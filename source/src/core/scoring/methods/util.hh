// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/methods/util.hh
/// @brief utility methods for scoring.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_methods_util_hh
#define INCLUDED_core_scoring_methods_util_hh

#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/pose/Pose.fwd.hh>

namespace core {
namespace scoring {
namespace methods {

core::Real get_residue_weight_by_ss(
	char ss
);

bool residues_interact(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::Real const interaction_cutoff
);

bool atoms_interact(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::id::AtomID const & id1,
	core::id::AtomID const & id2,
	core::Real const interaction_cutoff
);

/// @brief Given two residues that may or may not be connected, determine which of the two, if any,
/// is the lower one and which is the upper.
/// @details Inputs are rsd1 and rsd2; outputs are rsd1_is_lo and rsd2_is_lo.  Both will be false if
/// the residues aren't conventionally connected (i.e. the C of one connected to the N of the other).
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void determine_lo_and_hi_residues(
	core::pose::Pose const &pose,
	core::Size const rsd1,
	core::Size const rsd2,
	bool &res1_is_lo,
	bool &res2_is_lo
);

/// @brief Determines whether a long-range energies container exists in the pose energies object.  If not,
/// creates a new one and appends the score type to it, if necessary.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void create_long_range_energy_container(
	core::pose::Pose &pose,
	core::scoring::ScoreType const scoretype,
	core::scoring::methods::LongRangeEnergyType const lr_type
);

} // namespace methods
} // namespace scoring
} // namespace core

#endif
