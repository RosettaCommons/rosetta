// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_trials.hh
/// @brief  rotamer trials module header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_rotamer_trials_hh
#define INCLUDED_core_pack_rotamer_trials_hh

// pack headers
#include <core/pack/task/PackerTask.fwd.hh>

// pose headers
#include <core/pose/Pose.fwd.hh>

// scoring headers
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/types.hh>

// utility headers
#include <utility/vector1.hh>

namespace core {
namespace pack {

utility::vector1< uint >
repackable_residues( task::PackerTask const & the_task );


void
rotamer_trials(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP input_task
);

void
symmetric_rotamer_trials(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP input_task
);

}
}

#endif
