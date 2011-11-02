// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RationalMonteCarlo.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/moves/RationalMonteCarlo.hh>

// C/C++ headers
#include <iostream>
#include <string>

// External headers
#include <boost/function.hpp>
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/Tracer.hh>

// Project headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <protocols/viewer/viewers.hh>

// Package headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

using core::Real;
using core::Size;
using core::scoring::ScoreFunctionOP;
using protocols::moves::MoverOP;

typedef protocols::moves::Mover Parent;

static basic::Tracer TR("protocols.moves.RationalMonteCarlo");

RationalMonteCarlo::RationalMonteCarlo(MoverOP mover, ScoreFunctionOP score, Size num_trials, Real temperature, bool recover_low)
    : Parent("RationalMonteCarlo"), mover_(mover), num_trials_(num_trials), recover_low_(recover_low), next_trigger_id_(0) {
  mc_ = new protocols::moves::MonteCarlo(*score, temperature);
  protocols::viewer::add_monte_carlo_viewer(*mc_, "RationalMonteCarlo");
}

Size RationalMonteCarlo::num_trials() const {
  return num_trials_;
}

void RationalMonteCarlo::apply(core::pose::Pose& pose) {
  using core::pose::Pose;

  // Initialize the MonteCarlo object
  mc_->reset(pose);
  mc_->reset_counters();

  for (Size i = 1; i <= num_trials(); ++i) {
    // retain a copy of the pose in the event that the move is rejected
    Pose copy(pose);
    mover_->apply(pose);

    if (mc_->boltzmann(pose)) {  // accept
      fire_all_triggers(pose);
    } else {                     // reject
      pose = copy;
    }
  }

  // optionally recover the low-scoring pose
  if (recover_low())
    mc_->recover_low(pose);

  // show simulation statistics
  mc_->show_counters();
  mc_->score_function().show(TR, pose);
  TR.flush();
}

std::string RationalMonteCarlo::get_name() const {
  return "RationalMonteCarlo";
}

bool RationalMonteCarlo::recover_low() const {
  return recover_low_;
}

/// @detail Adds the specified trigger, returning a unique trigger id
int RationalMonteCarlo::add_trigger(const RationalMonteCarloTrigger& trigger) {
  const int tid = ++next_trigger_id_;
  triggers_[tid] = trigger;
  return tid;
}

/// @detail Attempts to remove the trigger with the given id. Issues a
/// warning if one is not found.
void RationalMonteCarlo::remove_trigger(int trigger_id) {
  Triggers::iterator i = triggers_.find(trigger_id);
  if (i == triggers_.end()) {
    TR.Warning << "Attempt to remove invalid trigger_id => " << trigger_id << std::endl;
    return;
  }
  triggers_.erase(i);
}

/// @detail Invokes all triggers registered with the pose
void RationalMonteCarlo::fire_all_triggers(const Pose& pose) {
  for (Triggers::iterator i = triggers_.begin(); i != triggers_.end(); ++i) {
    RationalMonteCarloTrigger& t = i->second;
    t(pose);
  }
}

}  // namespace moves
}  // namespace protocols
