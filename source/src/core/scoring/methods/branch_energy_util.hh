// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/branch_energyutil.hh
/// @brief  Utility functions for scoring branches.
/// @author Andrew Watkins

#ifndef INCLUDED_core_scoring_methods_branch_energy_util_hh
#define INCLUDED_core_scoring_methods_branch_energy_util_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>


namespace core {
namespace scoring {
namespace methods {


struct ResidueAtomOverlaps {

	Size res1;
	Size res2;

	std::string res1_ovl1_overlaps;
	std::string res1_ovl2_overlaps;

	std::string res2_ovu1_overlaps;
};

void
find_relevant_connections_onersd( pose::Pose const & pose, Size const seqpos, ResidueAtomOverlaps & branch_connection );

void
find_relevant_connections( pose::Pose const & pose, utility::vector1< ResidueAtomOverlaps > & branch_connections );

} // namespace methods
} // namespace scoring
} // namespace core

#endif
