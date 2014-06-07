// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/AbscriptLoopCloserCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_AbscriptLoopCloserCM_hh
#define INCLUDED_protocols_abinitio_abscript_AbscriptLoopCloserCM_hh

// Unit Headers
#include <protocols/abinitio/abscript/AbscriptLoopCloserCM.fwd.hh>
#include <protocols/environment/ClaimingMover.hh>

// Package headers

// Project headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/conformation/Conformation.fwd.hh>

#include <core/fragment/FragSet.hh>
#include <core/scoring/ScoreFunction.hh>

// C++ Headers

// ObjexxFCL Headers

//Req'd on WIN32
#include <basic/datacache/WriteableCacheableMap.hh>

namespace protocols {
namespace abinitio {
namespace abscript {

class AbscriptLoopCloserCM : public protocols::environment::ClaimingMover {
  typedef ClaimingMover Parent;
  typedef environment::claims::EnvClaims EnvClaims;

public:
  AbscriptLoopCloserCM();

  AbscriptLoopCloserCM( core::fragment::FragSetCOP fragset,
                        core::scoring::ScoreFunctionOP scorefxn );

  virtual ~AbscriptLoopCloserCM() {};

  virtual EnvClaims yield_claims( core::pose::Pose const&,
                                  basic::datacache::WriteableCacheableMapOP );

  virtual void broking_finished( environment::EnvClaimBroker::BrokerResult const& );

  std::string const& label() const { return label_; }

  void set_label( std::string const& label ) { label_ = label; }

  virtual std::string get_name() const;

  virtual void apply( core::pose::Pose& );

  virtual void
  parse_my_tag(utility::tag::TagCOP tag,
               basic::datacache::DataMap & data,
               protocols::filters::Filters_map const & filters,
               protocols::moves::Movers_map const & movers,
               core::pose::Pose const & pose );

  virtual
  moves::MoverOP clone() const;

protected:
  virtual void passport_updated();

private:
  void attempt_idealize( core::pose::Pose& );

  bool attempt_ccd( core::pose::Pose& );

  core::kinematics::FoldTreeOP final_ft_;
  core::fragment::FragSetCOP fragset_;
  mutable core::kinematics::MoveMapOP movemap_;
  core::scoring::ScoreFunctionOP scorefxn_;

  std::string label_;

  mutable bool bUpdateMM_;

}; // end AbscriptLoopCloserCM base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_AbscriptLoopCloserCM_hh
