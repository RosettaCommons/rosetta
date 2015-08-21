// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/seeded_abinitio/SeededAbinitio_util.cc
/// @brief
/// @author Eva-Maria Strauch ( evas01@u.washington.edu )


#ifndef INCLUDED_protocols_seeded_abinitio_SeededAbinitio_util_hh
#define INCLUDED_protocols_seeded_abinitio_SeededAbinitio_util_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
//#include <core/pack/task/TaskFactory.fwd.hh>
//#include <core/pack/task/operation/TaskOperation.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
//#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loops_main.hh>
//#include <protocols/loops/loops_main.hh>
//#include <protocols/loops/util.hh>
/////#include <protocols/loops/Loops.hh>
//#include <protocols/loops/LoopMover.fwd.hh>
//#include <protocols/loops/LoopMover.hh>


// Utillity Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility Headers
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace seeded_abinitio {

/// for runtime parsing of seeds
protocols::loops::Loops parse_seeds( core::pose::Pose const & pose, utility::vector1 < std::pair < std::string, std::string > > seed_vector);

/// to readjust the move map after length changes
void adjust_mm_to_length( core::pose::Pose const & pose, core::kinematics::MoveMapOP & mm );

void combine_two_poses( core::pose::Pose design_pose , core::pose::PoseOP target_chain );


}
}

#endif

