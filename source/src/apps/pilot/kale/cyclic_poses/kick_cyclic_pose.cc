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
using namespace protocols::loops::loop_closure;

// Global Variables {{{1

static thread_local basic::Tracer TR( "apps.pilot.kale" );

// Utility Functions {{{1

void describe_move(
		pose::Pose& model, 
		scoring::ScoreFunctionOP score_function) {

	using namespace std;
	int residues = model.total_residue();
	int width = 7 + 11 * residues - 1;

	TR << endl << string(width, '>');
	TR << endl << "score  " << setw(10) << score_function->score(model);

	TR << endl << "phi    ";
	for (int i = 1; i <= model.total_residue(); ++i) {
		TR << setw(10) << model.phi(i) << " ";
	}

	TR << endl << "psi    ";
	for (int i = 1; i <= model.total_residue(); ++i) {
		TR << setw(10) << model.psi(i) << " ";
	}

	TR << endl << "omega  ";
	for (int i = 1; i <= model.total_residue(); ++i) {
		TR << setw(10) << model.omega(i) << " ";
	}

	TR << endl << string(width, '<') << endl;
}

// Application {{{1

int main(int argc, char* argv []) {

		devel::init(argc, argv);

		// Define the variables used in this function.

		//Size first_index = 12, last_index = 2, cut_index = 14;
		Size first_index = 6, last_index = 2, cut_index = 4;
		std::string input_path = "structures/ideal_loop.pdb";
		std::string output_path = "structures/kicked_loop.pdb";

		pose::Pose model;
		scoring::ScoreFunctionOP score_function;
		kinematic_closure::KinematicMoverOP mover;
		kinematic_closure::TorsionSamplingKinematicPerturberOP perturber;

		// Initialize the variables used in this function.

		score_function = scoring::get_score_function();
		import_pose::pose_from_pdb(model, input_path);

		mover = new kinematic_closure::KinematicMover();
		perturber = new
			kinematic_closure::TorsionSamplingKinematicPerturber(mover.get());

		mover->set_pivots(first_index, cut_index, last_index);
		mover->set_perturber(perturber);

		// Attempt to make a kinematic closure move.

		mover->apply(model);

		// Write the resulting pose to a file.

		model.dump_pdb(output_path);
}

// }}}1
