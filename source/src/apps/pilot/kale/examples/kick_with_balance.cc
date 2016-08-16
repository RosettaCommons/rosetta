// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1

#include <devel/init.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// Namespaces {{{1
using namespace core;

#include <devel/balanced_kic/KinematicMover.hh>
using namespace devel::balanced_kic;

// #include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
// using namespace protocols::loops::loop_closure::kinematic_closure;

// Global Variables {{{1

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kale" );

// }}}1

int main(int argc, char* argv []) {

	devel::init(argc, argv);

	// Define the variables used in this function.

	pose::Pose pose;
	KinematicMoverOP mover;
	Size first_index, last_index, cut_index;
	Real temperature = 300;

	// Initialize the variables used in this function.

	import_pose::pose_from_file(pose, "structures/ideal_chain.7.pdb", core::import_pose::PDB_file);

	first_index = 2;
	cut_index = 4;
	last_index = 6;

	mover = new KinematicMover();
	mover->set_pivots(first_index, cut_index, last_index);
	mover->set_temperature(temperature);

	// Apply a kinematic closure move to the pose.  Because the move seems to be 
	// random (i.e. I can't fix the random seed), I can't compare future results 
	// against a single fixed structure.  However, I do know that this structure 
	// with these pivots should generally give good results.  Hopefully this will 
	// change if I accidentally break the closure algorithm.

	mover->apply(pose);

	TR << "move accepted? ";
	TR << (mover->last_move_succeeded() ? "yes" : "no") << std::endl;

	// Save the resulting structure.
	
	pose.dump_pdb("trial-closure.pdb");

}

