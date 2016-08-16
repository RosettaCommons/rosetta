// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/optimizeH.hh
/// @brief  declaration of standard hydrogen optimization subroutine
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_optimizeH_hh
#define INCLUDED_core_pack_optimizeH_hh

// Project headers
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {

/// @brief Call optimizeH and tell the user what chi angles have changed
void
optimize_H_and_notify(
	pose::Pose & pose,
	id::AtomID_Mask const & missing
);


/// @brief This function will optimize the placement of all movable
/// hydrogen atoms.  This includes the hydroxyl hydrogens as well as
/// the HIS protonation state.  If the -flip_HNQ flag is on the command
/// line, then it will also consider the flip states of histadine,
/// asparagine and glutamine, (nearly) as described by Word et al. 1999.
void
optimizeH(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn
);

}
}

#endif
