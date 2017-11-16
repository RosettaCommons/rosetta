// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/rational_mc/RationalMonteCarlo.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>

// C/C++ headers
#include <iostream>
#include <string>

// External headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/Tracer.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <protocols/viewer/viewers.hh>

// Package headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace simple_moves {
namespace rational_mc {

using core::Real;
using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using protocols::moves::Mover;
using protocols::moves::MoverOP;

static basic::Tracer TR( "protocols.simple_moves.RationalMonteCarlo" );

RationalMonteCarlo::RationalMonteCarlo(MoverOP mover, ScoreFunctionOP score, Size num_trials, Real temperature, bool recover_low)
: Mover("RationalMonteCarlo"), mover_(mover), num_trials_(num_trials), recover_low_(recover_low), next_trigger_id_(0) {
	mc_ = moves::MonteCarloOP( new protocols::moves::MonteCarlo(*score, temperature) );
	protocols::viewer::add_monte_carlo_viewer(*mc_, "RationalMonteCarlo");
}

void RationalMonteCarlo::apply(Pose& pose) {
	mc_->reset(pose);
	mc_->reset_counters();

	for ( Size i = 1; i <= num_trials(); ++i ) {
		Pose copy(pose);
		mover_->apply(pose);

		if ( mc_->boltzmann(pose) ) {  // accept
			fire_all_triggers(pose);
		} else {                     // reject
			pose = copy;
		}
	}

	if ( recover_low() ) {
		mc_->recover_low(pose);
	}

	// Show simulation statistics
	mc_->show_counters();
	mc_->score_function().show(TR, pose);
	TR.flush();
}

std::string RationalMonteCarlo::get_name() const {
	return "RationalMonteCarlo";
}

Size RationalMonteCarlo::add_trigger(const RationalMonteCarloTrigger& trigger) {
	Size tid = ++next_trigger_id_;
	triggers_[tid] = trigger;
	return tid;
}

void RationalMonteCarlo::remove_trigger(Size trigger_id) {
	Triggers::const_iterator i = triggers_.find(trigger_id);
	if ( i == triggers_.end() ) {
		TR.Warning << "Attempt to remove invalid trigger_id => " << trigger_id << std::endl;
		return;
	}
	triggers_.erase(i);
}

void RationalMonteCarlo::fire_all_triggers(const Pose& pose) {
	for ( Triggers::iterator i = triggers_.begin(); i != triggers_.end(); ++i ) {
		RationalMonteCarloTrigger& t = i->second;
		t(pose);
	}
}

// -- Accessors -- //
Size RationalMonteCarlo::num_trials() const {
	return num_trials_;
}

bool RationalMonteCarlo::recover_low() const {
	return recover_low_;
}

MoverOP RationalMonteCarlo::mover() const {
	return mover_;
}

Real RationalMonteCarlo::temperature() const {
	return mc_->temperature();
}

const ScoreFunction& RationalMonteCarlo::score_function() const {
	return mc_->score_function();
}

const Pose& RationalMonteCarlo::lowest_score_pose() const {
	return mc_->lowest_score_pose();
}

const Pose& RationalMonteCarlo::last_accepted_pose() const {
	return mc_->last_accepted_pose();
}

// -- Mutators -- //
void RationalMonteCarlo::set_num_trials(Size num_trials) {
	num_trials_ = num_trials;
}

void RationalMonteCarlo::set_recover_low(bool recover_low) {
	recover_low_ = recover_low;
}

void RationalMonteCarlo::set_mover(MoverOP mover) {
	mover_ = mover;
}

void RationalMonteCarlo::set_temperature(Real temperature) {
	mc_->set_temperature(temperature);
}

void RationalMonteCarlo::set_score_function(ScoreFunctionOP score) {
	mc_->score_function(*score);
}

void RationalMonteCarlo::reset(const Pose& pose) {
	mc_->reset(pose);
}

void RationalMonteCarlo::enable_autotemp(core::Real quench) {
	mc_->set_autotemp(true, quench);
}

void RationalMonteCarlo::disable_autotemp() {
	mc_->set_autotemp(false, 0);
}

}  // namespace rational_mc
}  // namespace simple_moves
}  // namespace protocols
