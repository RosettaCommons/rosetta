// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/inverse/util.hh
/// @brief  Utility functions for calculating jumps by knowing desired atom positions
/// @author Jack Maguire


#ifndef INCLUDED_core_kinematics_inverse_util_HH
#define INCLUDED_core_kinematics_inverse_util_HH


// Package headers
#include <core/kinematics/inverse/AlignmentAtom.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.fwd.hh>

namespace core {
namespace kinematics {
namespace inverse {

struct AlignmentAtomArray;

void
assert_atoms_are_downstream_of_jump(
	conformation::Conformation const & conformation,
	core::Size const jump_id,
	AlignmentAtomArray const & atom_arr
);

void
assert_atoms_are_upstream_of_jump(
	conformation::Conformation const & conformation,
	core::Size const jump_id,
	AlignmentAtomArray const & atom_arr
);

} // namespace inverse
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_inverse_util_HH
