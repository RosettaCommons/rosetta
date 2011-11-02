// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/min_pack.hh
/// @brief  pack and minimize sidechains
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_min_pack_hh
#define INCLUDED_core_pack_min_pack_hh

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {

/// @brief Interface function to the minimize-each-rotamer-during-packing sidechain placement algorithm.
void
min_pack(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task
);

void
stochastic_pack(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task
);

} // namespace pack
} // namespace core

#endif
