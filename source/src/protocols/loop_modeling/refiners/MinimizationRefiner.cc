// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefinerCreator.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

// Core headers
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace loop_modeling {
namespace refiners {

using namespace std;
using core::optimization::MinimizerOptions;
using core::optimization::MinimizerOptionsOP;
using core::optimization::MinimizerOptionsCOP;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::simple_moves::MinMover;
using protocols::simple_moves::MinMoverOP;
using protocols::simple_moves::symmetry::SymMinMover;

protocols::moves::MoverOP MinimizationRefinerCreator::create_mover() const {
	return protocols::moves::MoverOP( new MinimizationRefiner );
}

std::string MinimizationRefinerCreator::keyname() const {
	return "MinimizationRefiner";
}

MinimizationRefiner::MinimizationRefiner(
	bool cartesian, MinimizerOptionsOP options) {

	use_cartesian(cartesian);
	set_min_options(options);
}

void MinimizationRefiner::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose) {

	LoopMover::parse_my_tag(tag, data, filters, movers, pose);
	utilities::set_scorefxn_from_tag(*this, tag, data);
}

bool MinimizationRefiner::do_apply(Pose & pose) {
	using core::kinematics::MoveMap;
	using core::kinematics::MoveMapOP;
	using core::pose::symmetry::is_symmetric;
	using core::scoring::ScoreFunctionCOP;
	using protocols::loops::loops_set_move_map;

	pose.update_residue_neighbors();

	Loops const & loops = *get_loops();
	MinMoverOP minimizer( is_symmetric(pose) ? new SymMinMover : new MinMover );
	ScoreFunctionCOP score_function = get_score_function();
	MoveMapOP move_map( new MoveMap ); loops_set_move_map(
		pose, loops,
		/*fix sidechains:*/ false,
		*move_map,
		/*neighbor radius:*/ 10.0,
		/*allow omega moves:*/ true,
		/*allow takeoff torsion moves:*/ false);

	minimizer->score_function(score_function);
	minimizer->movemap(move_map);
	minimizer->min_options(min_options_);
	minimizer->cartesian(use_cartesian_);
	minimizer->apply(pose);

	return true;
}

ScoreFunctionOP MinimizationRefiner::get_score_function() {
	return get_tool<ScoreFunctionOP>(ToolboxKeys::SCOREFXN);
}

void MinimizationRefiner::set_score_function(ScoreFunctionOP score_function) {
	set_tool(ToolboxKeys::SCOREFXN, score_function);
}

void MinimizationRefiner::set_min_options(MinimizerOptionsOP options) {
	// If no minimizer options are given, use default values that seem to work
	// well.  These values were chosen based a simple benchmark run, but in that
	// benchmark no parameters really performed that much better or worse than
	// any others.  So loop modeling doesn't seem to be that sensitive to the
	// exact choice of minimizer options.

	if ( ! options ) {
		min_options_ = MinimizerOptionsOP( new MinimizerOptions(
			"lbfgs_armijo_nonmonotone",   // min_type
			1e-3,       // min_tolerance
			true,       // use_nblist
			false) );   // deriv_check
	} else {
		min_options_ = options;
	}
}

MinimizerOptionsOP MinimizationRefiner::get_min_options() {
	return min_options_;
}

MinimizerOptionsCOP MinimizationRefiner::get_min_options() const {
	return min_options_;
}

void MinimizationRefiner::use_cartesian(bool setting) {
	use_cartesian_ = setting;
}

bool MinimizationRefiner::use_cartesian() const {
	return use_cartesian_;
}

}
}
}
