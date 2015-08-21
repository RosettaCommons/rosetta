// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_indel_IndelOptimizationMover_fwd_hh
#define INCLUDED_protocols_indel_IndelOptimizationMover_fwd_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/indel/IndelOptimizationMover.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace indel {

class IndelOptimizationMover;

typedef utility::pointer::shared_ptr< IndelOptimizationMover > IndelOptimizationMoverOP;
typedef utility::pointer::shared_ptr< IndelOptimizationMover const > IndelOptimizationMoverCOP;


}
}

#endif
