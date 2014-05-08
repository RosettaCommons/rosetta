// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/FragmentCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/FragmentCM.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/CutBiasClaim.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>

#include <core/fragment/FragSet.hh>

#include <core/pose/Pose.hh>

#include <protocols/simple_moves/FragmentMover.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.FragmentCM", basic::t_info);

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace environment;

FragmentCM::FragmentCM():
  ClaimingMover(),
  mover_( NULL ),
  label_( "BASE" ),
  bUpdateMM_( true )
{}

FragmentCM::FragmentCM( simple_moves::FragmentMoverOP mover,
                        std::string const& label ):
  ClaimingMover(),
  label_( label )
{
  set_mover( mover );
}

FragmentCM::~FragmentCM() {}

void FragmentCM::set_label( std::string const& label ){
  label_ = label;
}

void FragmentCM::set_mover( simple_moves::FragmentMoverOP mover ){
  mover_ = mover;
  type( get_name() );
}


claims::EnvClaims FragmentCM::yield_claims( core::pose::Pose const&,
                                            basic::datacache::WriteableCacheableMapOP ){
  using namespace claims;
  claims::EnvClaims claim_list;

  TorsionClaimOP new_claim = new TorsionClaim( this, label(),
                                               std::make_pair( mover()->fragments()->min_pos(),
                                                               mover()->fragments()->max_pos() ) );

  new_claim->init_strength( DOES_NOT_INITIALIZE );
  new_claim->ctrl_strength( CAN_CONTROL );
  claim_list.push_back( new_claim );

  return claim_list;
}

void FragmentCM::initialize( Pose& pose ){
  if( bUpdateMM_ ){ update_movemap( pose ); }

  DofUnlock activation( pose.conformation(), passport() );
  mover()->apply_at_all_positions( pose );
}

void FragmentCM::apply( Pose& pose ) {
  if( bUpdateMM_ ){ update_movemap( pose ); }

  if( pose.conformation().is_protected() ) {
    DofUnlock activation( pose.conformation(), passport() );
    mover()->apply( pose );
  } else {
    mover()->apply( pose );
  }
}

void FragmentCM::update_movemap( Pose const& pose ) const {
  assert( mover() );
  if( has_passport() ){
    mover()->set_movemap( passport()->render_movemap( pose.conformation() ) );
  } else {
    core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap();
    mm->set_bb( false );
    mover()->set_movemap( mm );
  }
  bUpdateMM_ = false;
}

void FragmentCM::passport_updated() {
  bUpdateMM_ = true;
}

std::string FragmentCM::get_name() const {
  return "FragmentCM("+mover()->get_name()+")";
}


} // abscript
} // abinitio
} // protocols
