// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/LoopBuilder.hh>

// Core headers
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol headers
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/perturbers/IdealizeNonPhiPsi.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/EndToEndPivots.hh>
#include <protocols/kinematic_closure/pivot_pickers/LoopPivots.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>
#include <protocols/loop_modeling/refiners/LocalMinimizationRefiner.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <basic/Tracer.hh>

// Debug headers
#include <core/chemical/VariantType.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

namespace protocols {
namespace loop_modeling {

using namespace std;
using core::scoring::ScoreFunctionOP;

LoopBuilder::LoopBuilder() { // {{{1
	using protocols::kinematic_closure::KicMover;
	using protocols::kinematic_closure::perturbers::IdealizeNonPhiPsi;
	using protocols::kinematic_closure::perturbers::RamaPerturber;
	using protocols::kinematic_closure::pivot_pickers::EndToEndPivots;
	using protocols::kinematic_closure::pivot_pickers::LoopPivots;
	using protocols::kinematic_closure::solution_pickers::FilteredSolutions;
	using protocols::kinematic_closure::solution_pickers::FilteredSolutionsOP;
	using protocols::loop_modeling::refiners::LocalMinimizationRefiner;

	// The rama check works by comparing the generated torsions to the input
	// torsions.  Since the purpose of loop rebuilding is to forget everything
	// about the input structure, the rama check shouldn't be used.

	FilteredSolutionsOP solution_picker = new FilteredSolutions;
	solution_picker->dont_check_rama();
	solution_picker->be_lenient();

	kic_mover_ = new KicMover;
	kic_mover_->add_perturber(new IdealizeNonPhiPsi);
	kic_mover_->add_perturber(new RamaPerturber);
	kic_mover_->set_pivot_picker(new LoopPivots);
	kic_mover_->set_solution_picker(solution_picker);

	minimizer_ = new LocalMinimizationRefiner;

	register_nested_loop_mover(kic_mover_);
	register_nested_loop_mover(minimizer_);
}

bool LoopBuilder::do_apply(Pose & pose, Loop const & loop) { // {{{1
	basic::Tracer tr("protocols.loop_modeling.LoopBuilder");

	// Only attempt to rebuild loops that are marked as "extended".

	if (! loop.is_extended()) return true;

	// Setup the loop movers.

	kic_mover_->set_loop(loop);
	minimizer_->set_loop(loop);

	// Make a strong effort to rebuild the loop with KIC.

	for (Size i = 1; i <= 1000 && ! kic_mover_->was_successful(); i++) {
		tr << "Loop building attempt: " << i << endl;
		kic_mover_->apply(pose);
	}

	if (! kic_mover_->was_successful()) return false;

	// Minimize the loop if it was successfully built.

	minimizer_->apply(pose);
	return minimizer_->was_successful();
}
// }}}1

}
}
