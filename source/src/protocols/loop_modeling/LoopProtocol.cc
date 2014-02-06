// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/LoopBuilder.hh>
#include <protocols/loop_modeling/LoopMover.hh>

// Protocols headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility headers
#include <utility/exit.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <cmath>

#define foreach BOOST_FOREACH
using namespace std;

namespace protocols {
namespace loop_modeling {

using core::Real;
using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::loops::Loop;
using protocols::loops::Loops;

LoopProtocol::LoopProtocol() { // {{{1
	mover_ = NULL;
	monte_carlo_ = NULL;
	score_function_ = NULL;

	iterations_ = IndexList(3, 1);
	ramp_score_function_ = false;
	ramp_temperature_ = true;

	initial_temperature_ = 1.5;
	final_temperature_ = 0.5;
}

LoopProtocol::~LoopProtocol() {} // {{{1

// }}}1

void LoopProtocol::apply(Pose & pose) { // {{{1
	start_protocol(pose);

	for (Size i = 1; i <= iterations_[1]; i++) {
		ramp_score_function(i);

		for (Size j = 1; j <= iterations_[2]; j++) {
			ramp_temperature(j);

			for (Size k = 1; k <= iterations_[3]; k++) {
				attempt_loop_move(pose, i, j, k);
			}
		}
	}

	finish_protocol(pose);
}

void LoopProtocol::start_protocol(Pose & pose) { // {{{1
	if (! loop_.start() || ! loop_.stop()) {
		utility_exit_with_message("No loop region specified.");
	}
	if (! mover_) {
		utility_exit_with_message("No LoopMover specified.");
	}
	if (! score_function_) {
		score_function_ = core::scoring::getScoreFunction();
	}

	Pose native = pose;
	pose.update_residue_neighbors();

	// Rebuild the loop from scratch to avoid bias.
	//LoopBuilder builder(loop_, score_function_);
	//builder.apply(pose);

	// Setup the loop mover.  This must be done before the Monte Carlo object is 
	// initialized, because otherwise the "last accepted pose" in the Monte Carlo 
	// object will not have been properly setup.  In particular, it will have the 
	// wrong fold tree. 
	
	mover_->set_loop(loop_);
	mover_->set_score_function(score_function_);
	mover_->setup(pose);

	// Setup the loggers.
	foreach (loggers::LoggerOP logger, loggers_) {
		logger->log_beginning(pose, iterations_);
		mover_->add_logger(logger);
	}

	// Setup the acceptance criterion and the temperature.
	monte_carlo_ = new protocols::moves::MonteCarlo(
			pose, *score_function_, initial_temperature_);

	temperature_scale_factor_ = std::pow(
			final_temperature_ / initial_temperature_,
			1.0 / (iterations_[1] * iterations_[2]));

	// Report the initial score and RMSD.
	//using protocols::loops::loop_rmsd;
	//Loops loop_wrapper;
	//loop_wrapper.add_loop(loop_);

	//cout << "Initial RMSD:  " << loop_rmsd(native, pose, loop_wrapper) << endl;
	//cout << "Initial Score: " << score_function_->score(pose) << endl;
	//cout << endl;
}

void LoopProtocol::ramp_score_function(Size iteration) { // {{{1
	using core::scoring::fa_rep;
	using core::scoring::rama;

	if (ramp_score_function_) {
		Real ramp_factor = 1 / (iteration - iterations_[1] + 1);

		Real repulsive_weight = score_function_->get_weight(fa_rep) * ramp_factor;
		Real rama_weight = score_function_->get_weight(rama) * ramp_factor;

		score_function_->set_weight(fa_rep, repulsive_weight);
		score_function_->set_weight(rama, rama_weight);
	}
}

void LoopProtocol::ramp_temperature(Size /*iteration*/) { // {{{1
	if (ramp_temperature_) {
		Real temperature = monte_carlo_->temperature();
		monte_carlo_->set_temperature(temperature * temperature_scale_factor_);
	}
}

void LoopProtocol::attempt_loop_move( // {{{1
		Pose & pose, Size i, Size j, Size k) {

	foreach (loggers::LoggerOP logger, loggers_) {
		 logger->log_iteration(pose, i, j, k);
	}

	mover_->apply(pose);
	monte_carlo_->boltzmann(pose);

	foreach (loggers::LoggerOP logger, loggers_) {
		 logger->log_monte_carlo(monte_carlo_);
	}
}

void LoopProtocol::finish_protocol(Pose const & pose) { // {{{1
	foreach (loggers::LoggerOP logger, loggers_) {
		 logger->log_ending(pose);
	}
}
// }}}1

void LoopProtocol::set_mover(LoopMoverOP mover) { // {{{1
	mover_ = mover;
}

void LoopProtocol::set_loop(Loop const & loop) { // {{{1
	loop_ = loop;
}

void LoopProtocol::set_score_function(ScoreFunctionOP function) { // {{{1
	score_function_ = function;
}

void LoopProtocol::set_iterations(IndexList iterations) { // {{{1
	iterations_ = iterations;
}

void LoopProtocol::set_iterations(Size i, Size j, Size k) { // {{{1
	iterations_.resize(3);
	iterations_[1] = i;
	iterations_[2] = j;
	iterations_[3] = k;
}

// This is a comment 

void LoopProtocol::set_temperature_schedule( // {{{1
		Real initial, Real final) {

	initial_temperature_ = initial;
	final_temperature_ = final;
}

void LoopProtocol::set_temperature_ramping(bool value) { // {{{1
	ramp_temperature_ = value;
}

void LoopProtocol::set_score_function_ramping(bool value) { // {{{1
	ramp_score_function_ = value;
}

void LoopProtocol::add_logger(loggers::LoggerOP logger) { // {{{1
	loggers_.push_back(logger);
}
// }}}1

}
}

