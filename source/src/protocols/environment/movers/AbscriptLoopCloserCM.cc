// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/AbscriptLoopCloserCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/movers/AbscriptLoopCloserCM.hh>

// Package headers
#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/EnvExcn.hh>

#include <core/environment/DofPassport.hh>

// Project headers
#include <protocols/loops/loop_closure/ccd/WidthFirstSlidingWindowLoopClosure.hh>
#include <protocols/loops/Exceptions.hh>

#include <protocols/checkpoint/CheckPointer.hh>

#include <protocols/jumping/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <protocols/idealize/IdealizeMover.hh>

//Utility Headers
#include <utility/excn/Exceptions.hh>

#include <boost/functional/hash.hpp>
#include <boost/foreach.hpp>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.AbscriptLoopCloserCM", basic::t_info);

namespace protocols {
namespace environment {

using namespace core::environment;

AbscriptLoopCloserCM::AbscriptLoopCloserCM( core::fragment::FragSetCOP fragset,
                                            core::scoring::ScoreFunctionOP scorefxn ):
  Parent(),
  fragset_( fragset ),
  scorefxn_( scorefxn )
{}

claims::EnvClaims AbscriptLoopCloserCM::yield_claims( core::pose::Pose& in_pose ){
  claims::EnvClaims claims;

  // We want to control everything that will be relevant to the output pose (which should be the same as the input.
  claims::TorsionClaimOP claim = new claims::TorsionClaim( this, "BASE", std::make_pair( 1, in_pose.total_residue() ) );
  claim->ctrl_strength( claims::CAN_CONTROL );
  claim->init_strength( claims::DOES_NOT_INITIALIZE );

  claims.push_back( claim );

  return claims;
}

void AbscriptLoopCloserCM::apply( core::pose::Pose& in_pose ){

  //Produce unprotected pose
  core::pose::Pose pose( in_pose );
  pose.set_new_conformation( new core::conformation::Conformation( in_pose.conformation() ) );

  tr.Debug << "Closing loops with current fold tree: " << in_pose.fold_tree();
  tr.Debug << "Using final fold tree: " << *final_ft_ << std::endl;

  // Close Loops
  bool success;

  attempt_idealize( pose );

  success = attempt_ccd( pose );

  if( success ){
    attempt_idealize( pose );
  } else {
    // We could implement "don't fail unclosed" here, but I don't think we need to?
    return;
  }

  // Copy result into protected conformation in in_pose
  DofUnlock unlock( in_pose.conformation(), passport() );

  for ( Size i = 1; i <= in_pose.total_residue(); ++i ) {
    try {
      if( pose.omega( i ) != in_pose.omega( i ) ){
        in_pose.set_omega( i, pose.omega( i ) );
      }
      if( pose.phi( i ) != in_pose.phi( i ) ){
        in_pose.set_phi( i, pose.phi( i ) );
      }
      if( pose.psi( i ) != in_pose.psi( i ) ) {
        in_pose.set_psi( i, pose.psi( i ) );
      }
    } catch ( EXCN_Env_Security_Exception& ){
      tr.Error << "[ERROR] Unauthorized changes occurred during loop closure (attempt to write to resid "
               << i << ")." << std::endl;
      throw;
    }
  }
}

bool AbscriptLoopCloserCM::attempt_ccd( core::pose::Pose& pose ){
  using namespace loops::loop_closure::ccd;

  try {
    checkpoint::CheckPointer checkpointer( "EnvAbscriptLoopCloserCM" );
    SlidingWindowLoopClosureOP closing_protocol;
    closing_protocol = new WidthFirstSlidingWindowLoopClosure( fragset_,
                                                               scorefxn_,
                                                               movemap_ );

    jumping::close_chainbreaks( closing_protocol,
                               pose,
                               checkpointer,
                               get_current_tag(),
                               *final_ft_ );
  } catch ( loops::EXCN_Loop_not_closed& excn ) {
    set_last_move_status( moves::FAIL_DO_NOT_RETRY );
    set_current_tag( "C_"+get_current_tag().substr(std::min(2,(int)get_current_tag().size())) );
    return false;
  }

  return true;
}

void AbscriptLoopCloserCM::attempt_idealize( core::pose::Pose& pose ) {
  idealize::IdealizeMoverOP idealizer(new idealize::IdealizeMover);
  idealizer->fast( false );

  utility::vector1< core::Size > pos_list;
  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    if( movemap_->get_bb( i ) ){
      pos_list.push_back( i );
    }
  }

  idealizer->set_pos_list( pos_list );

  core::kinematics::FoldTree fold_tree( pose.fold_tree() );
  pose.fold_tree( *final_ft_ );
  idealizer->apply( pose );
  pose.fold_tree( fold_tree );
  (*scorefxn_)( pose );
  //jd2::output_intermediate_pose( pose, "loops_closed_preprocessed" );
}

void AbscriptLoopCloserCM::passport_updated() {
  if( has_passport() ){
    movemap_ = passport()->render_movemap();
  } else {
    movemap_ = new core::kinematics::MoveMap();
    movemap_->set_bb( false );
  }
}


void AbscriptLoopCloserCM::broking_finished( EnvClaimBroker::BrokerResult const& broker ) {
  //set the target fold tree for the closer.
  final_ft_ = broker.closer_ft;
}

std::string AbscriptLoopCloserCM::get_name() const {
  return "AbscriptLoopCloserCM";
}

} // environment
} // protocols
