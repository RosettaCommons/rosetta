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
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.hh>

// Protocols headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
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
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::moves::MonteCarloOP;

LoopProtocol::LoopProtocol() { // {{{1
	movers_ = new utilities::LoopMoverGroup;
	register_nested_loop_mover(movers_);

	monte_carlo_ = NULL;

	iterations_ = IndexList(3, 1);
	ramp_score_function_ = false;
	ramp_temperature_ = true;

	initial_temperature_ = 1.5;
	final_temperature_ = 0.5;
}

LoopProtocol::~LoopProtocol() {} // {{{1

// }}}1

bool LoopProtocol::do_apply(Pose & pose) { // {{{1
	start_protocol(pose);

	for (Size i = 1; i <= iterations_[1]; i++) {
		ramp_score_function(i);
		monte_carlo_->recover_low(pose);

		for (Size j = 1; j <= iterations_[2]; j++) {
			ramp_temperature(j);

			for (Size k = 1; k <= iterations_[3]; k++) {
				attempt_loop_move(pose, i, j, k);
			}
		}
	}

	pose = monte_carlo_->lowest_score_pose();
	finish_protocol(pose);

	return true;
}

void LoopProtocol::start_protocol(Pose & pose) { // {{{1
	if (movers_->empty()) {
		utility_exit_with_message("No movers specified.");
	}

	// Setup the Monte Carlo simulation.

	protocols::loops::add_cutpoint_variants(pose);
	protocols::loops::loop_mover::loops_set_chainbreak_weight(
			get_score_function(), 1);

	monte_carlo_ = new protocols::moves::MonteCarlo(
			pose, *get_score_function(), initial_temperature_);

	temperature_scale_factor_ = std::pow(
			final_temperature_ / initial_temperature_,
			1.0 / (iterations_[1] * iterations_[2]));

	// Setup the loggers.

	foreach (loggers::LoggerOP logger, loggers_) {
		logger->log_beginning(pose, iterations_);
	}
}

void LoopProtocol::ramp_score_function(Size iteration) { // {{{1
	using core::scoring::fa_rep;
	using core::scoring::rama;

	if (ramp_score_function_) {
		Real ramp_factor = 1 / (iteration - iterations_[1] + 1);

		Real repulsive_weight = get_score_function()->get_weight(fa_rep) * ramp_factor;
		Real rama_weight = get_score_function()->get_weight(rama) * ramp_factor;

		get_score_function()->set_weight(fa_rep, repulsive_weight);
		get_score_function()->set_weight(rama, rama_weight);
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

	movers_->apply(pose);

	if (movers_->was_successful()) {
		monte_carlo_->boltzmann(pose);
	} else {
		monte_carlo_->reset_last_accepted(pose);
	}

	foreach (loggers::LoggerOP logger, loggers_) {
		 logger->log_monte_carlo(monte_carlo_);
	}
}

void LoopProtocol::finish_protocol(Pose & pose) { // {{{1

	// Clean up the loggers.

	foreach (loggers::LoggerOP logger, loggers_) {
		 logger->log_ending(pose);
	}
}
// }}}1

void LoopProtocol::add_mover(LoopMoverOP mover) { // {{{1
	movers_->add_mover(mover);
}

void LoopProtocol::add_filter(protocols::filters::FilterOP filter) { // {{{1
	movers_->add_filter(filter);
}

void LoopProtocol::add_acceptance_check(string name) { // {{{1
	movers_->add_mover(new utilities::AcceptanceCheck(monte_carlo_, name));
}


void LoopProtocol::add_logger(loggers::LoggerOP logger) { // {{{1
	loggers_.push_back(logger);
}
// }}}1

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
// }}}1

}
}

