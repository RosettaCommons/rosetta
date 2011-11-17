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
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <protocols/viewer/viewers.hh>

// Package headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace moves {

using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::moves::MoverOP;

typedef protocols::moves::Mover Parent;

static basic::Tracer TR("protocols.moves.RationalMonteCarlo");


RationalMonteCarlo::RationalMonteCarlo(MoverOP mover, ScoreFunctionOP score, unsigned num_trials, double temperature, bool recover_low)
    : Parent("RationalMonteCarlo"), mover_(mover), num_trials_(num_trials), recover_low_(recover_low), next_trigger_id_(0) {
  mc_ = new protocols::moves::MonteCarlo(*score, temperature);
  protocols::viewer::add_monte_carlo_viewer(*mc_, "RationalMonteCarlo");
}

unsigned RationalMonteCarlo::num_trials() const {
  return num_trials_;
}

void RationalMonteCarlo::apply(Pose& pose) {
  mc_->reset(pose);
  mc_->reset_counters();

  for (unsigned i = 1; i <= num_trials(); ++i) {
    Pose copy(pose);
    mover_->apply(pose);

    if (mc_->boltzmann(pose)) {  // accept
      fire_all_triggers(pose);
    } else {                     // reject
      pose = copy;
    }
  }

  if (recover_low())
    mc_->recover_low(pose);

  // Show simulation statistics
  mc_->show_counters();
  mc_->score_function().show(TR, pose);
  TR.flush();
}

void RationalMonteCarlo::set_native(const Pose& native) {
  native_ = native;
  gdtmms_.clear();
  rmsds_.clear();
  add_trigger(boost::bind(&protocols::moves::RationalMonteCarlo::compute_analytics, this, _1));
}

const utility::vector1<double>& RationalMonteCarlo::gdtmms() const {
  return gdtmms_;
}

const utility::vector1<double>& RationalMonteCarlo::rmsds() const {
  return rmsds_;
}

std::string RationalMonteCarlo::get_name() const {
  return "RationalMonteCarlo";
}

bool RationalMonteCarlo::recover_low() const {
  return recover_low_;
}

int RationalMonteCarlo::add_trigger(const RationalMonteCarloTrigger& trigger) {
  unsigned tid = ++next_trigger_id_;
  triggers_[tid] = trigger;
  return tid;
}

void RationalMonteCarlo::remove_trigger(int trigger_id) {
  Triggers::const_iterator i = triggers_.find(trigger_id);
  if (i == triggers_.end()) {
    TR.Warning << "Attempt to remove invalid trigger_id => " << trigger_id << std::endl;
    return;
  }
  triggers_.erase(i);
}

void RationalMonteCarlo::fire_all_triggers(const Pose& pose) {
  for (Triggers::iterator i = triggers_.begin(); i != triggers_.end(); ++i) {
    RationalMonteCarloTrigger& t = i->second;
    t(pose);
  }
}

void RationalMonteCarlo::compute_analytics(const Pose& pose) {
  rmsds_.push_back(core::scoring::native_CA_rmsd(native_, pose));
  gdtmms_.push_back(core::scoring::native_CA_gdtmm(native_, pose));
}

}  // namespace moves
}  // namespace protocols
