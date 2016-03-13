// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ProlineFixMover.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/simple_moves/ProlineFixMover.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

// tracer
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

void ProlineFixMover::apply( core::pose::Pose & pose ) {

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );

	// try a repack
	core::pack::task::PackerTaskOP task(
		core::pack::task::TaskFactory::create_packer_task( pose )
	);
	task->initialize_from_command_line();
	task->restrict_to_repacking();
	core::pack::pack_rotamers(pose, (*scorefxn), task);

	// try a minimize
	core::kinematics::MoveMap mm;
	mm.set_bb(false);
	mm.set_chi(true);

	core::optimization::AtomTreeMinimizer().run(
		pose, mm, (*scorefxn), core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true)
	);

} // apply

std::string
ProlineFixMover::get_name() const {
	return "ProlineFixMover";
}

} // moves
} // protocols

