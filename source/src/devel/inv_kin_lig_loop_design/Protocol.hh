// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Protocol.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_PROTOCOL_HH
#define DEVEL_INVKINLIGLOOPDESIGN_PROTOCOL_HH

#include <vector>

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>

#include <devel/inv_kin_lig_loop_design/Loop.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace devel {

namespace inv_kin_lig_loop_design {

using namespace std;
using platform::Size;

// ==================================================
// ==================== Protocol ====================
// ==================================================

struct Protocol {
	Protocol() {}

	void setParams( TagCOP tag, vector<Loop> const& loops_ );
	void apply(core::pose::PoseOP pose);

private:

	int max_cycles;
	int max_cycles_start;
	int max_cycles_start_3mer;
	int max_cycles_start_3mer_anchor;
	int max_cycles_start_3mer_1mer;
	int max_cycles_start_1mer;
	int max_cycles_ramp;
	int max_cycles_ramp_sm;
	int max_cycles_ramp_sm_min;
	int max_cycles_design;
	int max_cycles_design_sm;
	int max_cycles_design_sm_min;

	core::pose::PoseOP pose;
	core::kinematics::MoveMapOP move_map;
	core::scoring::ScoreFunctionOP score_fxn_lores;
	core::scoring::ScoreFunctionOP score_fxn_hires;
	vector<Loop> loops;

	void phase_lores();
	void phase_hires();

};


}

}

#endif // DEVEL_LOOPDESIGN_PROTOCOL_HH
