// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/AbscriptStageMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_AbscriptStageMover_hh
#define INCLUDED_protocols_environment_AbscriptStageMover_hh

// Unit Headers
#include <protocols/environment/movers/AbscriptStageMover.fwd.hh>

// Package headers
#include <protocols/environment/movers/StageID.hh>
#include <protocols/environment/movers/AbscriptMover.hh>
#include <protocols/environment/movers/StagePreparer.fwd.hh>

#include <protocols/environment/ClaimingMover.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

// Project headers
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>

// Utility Headers
#include <utility/vector0.fwd.hh>

// C++ Headers
#include <set>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class AbscriptStageMover : public moves::Mover {
  typedef std::set< ClaimingMoverOP > MoverSet;
  typedef std::set< StagePreparerOP > PreparerSet;

public:
  AbscriptStageMover( StageID stage,
                      moves::MonteCarloOP mc,
                      core::scoring::ScoreFunctionOP score,
                      core::Size cycles );

  void apply( core::pose::Pose& );

  void yield_submovers( MoverSet& set ) const;

  std::string get_name() const;

  //@returns if this stage step should be run or not (sometimes they can be skipped).
  bool setup_stage( core::pose::Pose& pose, core::Real const& progress );

  void add_submover( ClaimingMoverOP mover, core::Real weight );

  void add_preparer( StagePreparerOP mover );

  void set_cycles_adjust( core::Real in );

  void set_seq_sep_ramping( core::Real const slope, core::Real const intercept );

  void set_chainbreak_ramping( core::Real const slope, core::Real const intercept );

  core::Real const& seq_sep_intercept() { return seqsep_intcpt_; }

  core::scoring::ScoreFunctionCOP scorefxn() const { return score_; }

private:

  void update_max_seq_sep( core::pose::Pose& pose, core::Real const& progress );

  StageID stage_;
  moves::RandomMoverOP random_mover_;
  MoverSet movers_;
  PreparerSet preparers_;
  core::Size const cycles_;
  core::Real cycles_adjust_;
  core::scoring::ScoreFunctionOP score_;
  moves::MonteCarloOP mc_;

  core::Real temperature_;

  core::Real seqsep_slope_;
  core::Real seqsep_intcpt_;
  core::Real chbreak_slope_;
  core::Real chbreak_intcpt_;

  constraints_additional::MaxSeqSepConstraintSetOP constraints_;
};

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_AbscriptStageMover_hh
