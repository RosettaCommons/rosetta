// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/AbscriptLoopCloserCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_AbscriptLoopCloserCM_hh
#define INCLUDED_protocols_environment_AbscriptLoopCloserCM_hh

// Unit Headers
#include <protocols/environment/movers/AbscriptLoopCloserCM.fwd.hh>
#include <protocols/environment/ClaimingMover.hh>

// Package headers

// Project headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class AbscriptLoopCloserCM : public ClaimingMover {
  typedef ClaimingMover Parent;
public:
  AbscriptLoopCloserCM( core::fragment::FragSetCOP fragset,
                core::scoring::ScoreFunctionOP scorefxn );

  virtual ~AbscriptLoopCloserCM() {};

  virtual claims::EnvClaims yield_claims( core::pose::Pose& );

  virtual void broking_finished( EnvClaimBroker::BrokerResult const& );

  virtual void passport_updated();

  virtual std::string get_name() const;

  virtual void apply( core::pose::Pose& );

  virtual bool is_loop_closer() const { return true; }

private:
  void attempt_idealize( core::pose::Pose& );

  bool attempt_ccd( core::pose::Pose& );

  core::kinematics::FoldTreeOP final_ft_;
  core::fragment::FragSetCOP fragset_;
  core::kinematics::MoveMapOP movemap_;
  core::scoring::ScoreFunctionOP scorefxn_;

}; // end AbscriptLoopCloserCM base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_AbscriptLoopCloserCM_hh
