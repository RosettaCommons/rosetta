// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1

#include <devel/init.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

#include <iomanip>

// Namespaces {{{1

using namespace core;
using namespace basic::options;
using namespace protocols::loops::loop_closure;

// Global Variables {{{1

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kale" );

// Application {{{1

int main(int argc, char* argv []) {

	devel::init(argc, argv);

	pose::Pose model;
	std::string input_path, output_path;

	// Make sure the needed command line flags were given.

	if (option[OptionKeys::in::file::s].active() == false) {
			utility_exit_with_message("No input PDB file specified.");
	}

	if (option[OptionKeys::loops::output_pdb].active() == false) {
			utility_exit_with_message("No output PDB file specified.");
	}

	input_path = option[OptionKeys::in::file::s]()[1];
	output_path = option[OptionKeys::loops::output_pdb]().name();

	// Read the pose into memory and immediately write it out again.

	import_pose::pose_from_pdb(model, input_path);
	model.dump_pdb(output_path);
}

// }}}1
