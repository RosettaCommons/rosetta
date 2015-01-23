// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/claims/EnvClaim.fwd.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>

#include <protocols/environment/claims/BrokerElements.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvClaimBroker.hh>

//Other headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/types.hh>

#include <test/core/init_util.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>
#include <core/pose/datacache/CacheableDataType.hh>

//C++ headers
#include <iostream>

using namespace protocols::environment;

// ---------------- Toy Movers --------------- //

namespace protocols {
namespace environment {

const core::Size CLAIMED_RESID = 10;
const core::Size UNCLAIMED_RESID = 9;
const core::Real NEW_PHI = 23.0;

const core::Size JUMP_START = 3;
const core::Size JUMP_END = 7;
class ToyMover : public protocols::environment::ClientMover {
protected:

  ToyMover( bool claim, bool move ):
    claim_( claim ),
    move_( move )
  {}

  //Yes I know this is against the law
  bool claim_;
  bool move_;
};

class TorsionMover : public ToyMover {
public:
  TorsionMover( bool claim, bool move,
    claims::ControlStrength c_str = claims::MUST_CONTROL,
    claims::ControlStrength i_str = claims::DOES_NOT_CONTROL ):
    ToyMover( claim, move ),
    resnum_( CLAIMED_RESID ),
    control_str_( c_str ),
    init_str_( i_str )
  {}

	using ToyMover::apply;
  void apply( Pose& pose, Size resid ){
    DofUnlock activation( pose.conformation(), passport() );
    if( move_ ) {
      pose.set_phi( resid, NEW_PHI );
    }
  }

  virtual void apply( Pose& pose ){
    apply( pose, resnum_ );
  }

  virtual void missing_unlock_apply( Pose& pose ){
    if( move_ ){
      pose.set_phi( resnum_, NEW_PHI );
    }
  }

  virtual std::string get_name() const{
    return "TorsionMover";
  }

  virtual claims::EnvClaims yield_claims( core::pose::Pose const&,
            basic::datacache::WriteableCacheableMapOP ){
    using core::environment::LocalPosition;
    claims::EnvClaims claims;

    if( claim_ ){
      claims::TorsionClaimOP new_claim( new claims::TorsionClaim(
		utility::pointer::dynamic_pointer_cast< protocols::environment::ClientMover >(get_self_ptr()),
		LocalPosition( "BASE", resnum_ ) ) );
      new_claim->strength( control_str_, init_str_ );
      claims.push_back( new_claim );
    }

    return claims;
  }

  virtual void initialize( core::pose::Pose& pose ){
    DofUnlock activation( pose.conformation(), passport() );
    if( move_ ){
      pose.set_phi( resnum_, NEW_PHI );
    }
  }

  Size resnum() {
    return resnum_;
  }

private:
  Size resnum_;
  claims::ControlStrength control_str_;
  claims::ControlStrength init_str_;
};

typedef utility::pointer::shared_ptr< TorsionMover > TorsionMoverOP;

} //environment
} //protocols

// --------------- Test Class --------------- //

class EnvironmentTorsion : public CxxTest::TestSuite {
public:

  // Shared data elements go here.
  core::pose::Pose pose;
  utility::vector1< core::Real > init_phis;

  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
    core_init();

    using namespace protocols::environment;
    using namespace core::environment;

    core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

    //store phi values for comparison in the test
    for( core::Size i = 1; i <= pose.total_residue(); ++i ){
      init_phis.push_back( pose.phi( i ) );
    }

  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_empty_environment(){
    TS_TRACE( "Beginning: test_empty_environment" );

    using namespace protocols::environment;
    using namespace core::environment;

    TorsionMoverOP mover( new TorsionMover( false, false ) );
    EnvironmentOP env_op( new Environment( "env" ) );
    Environment & env = *env_op;

    core::pose::Pose brokered;
    TS_ASSERT_THROWS_NOTHING( brokered = env.start( pose ) );
    core::pose::Pose done;
    TS_ASSERT_THROWS_NOTHING( done = env.end( brokered ) );

    TS_TRACE( "End: test_empty_environment" );
  }

  void test_dual_environment(){
    TS_TRACE( "Beginning: test_dual_environment" );

    using namespace protocols::environment;
    using namespace core::environment;

    TorsionMoverOP tier1_mover( new TorsionMover( true, true ) );
    TorsionMoverOP tier2_mover( new TorsionMover( true, true ) );

    EnvironmentOP env1_op( new Environment( "tier1" ) );
    Environment & env1 = *env1_op;
    EnvironmentOP env2_op( new Environment( "tier2" ) );
    Environment & env2 = *env2_op;

    env1.register_mover( tier1_mover );
    env2.register_mover( tier2_mover );

    core::pose::Pose final_pose;

    {
      core::pose::Pose protected1;
      TS_ASSERT_THROWS_NOTHING( protected1 = env1.start( pose ) );
      TS_ASSERT_THROWS( tier2_mover->apply( protected1 ), utility::excn::EXCN_NullPointer );
      TS_ASSERT_EQUALS( protected1.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ] );
      {
  core::pose::Pose protected2;
  TS_ASSERT_THROWS_NOTHING( protected2 = env2.start( protected1 ) );
  TS_ASSERT_THROWS( tier1_mover->apply( protected2 ), EXCN_Env_Security_Exception );
  TS_ASSERT_DELTA( protected2.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 0.001 );
  TS_ASSERT_THROWS_NOTHING( tier2_mover->apply( protected2 ) );
  TS_ASSERT_DELTA( protected2.phi( CLAIMED_RESID ), NEW_PHI, 0.001 );
  TS_ASSERT_THROWS_NOTHING( protected1 = env2.end( protected2 ) );
      }
      TS_ASSERT_DELTA( protected1.phi( CLAIMED_RESID ), NEW_PHI, 0.001 );
      TS_ASSERT_THROWS_NOTHING( tier1_mover->apply( protected1) );
      TS_ASSERT_THROWS( tier2_mover->apply( protected1 ), utility::excn::EXCN_NullPointer );
      TS_ASSERT_THROWS_NOTHING( final_pose = env1.end( protected1 ) );
    }

    for( core::Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
      if( seqpos == CLAIMED_RESID ){
  TS_ASSERT_DELTA( final_pose.phi( seqpos ), NEW_PHI, 0.000001 );
      } else {
  TS_ASSERT_DELTA( final_pose.phi( seqpos ), init_phis[ seqpos ], 0.000001 );
      }
    }

    TS_TRACE( "End: test_dual_environment" );
  }

  void test_single_phi_moves(){
    TS_TRACE( "Beginning: test_single_phi_moves" );
    using namespace protocols::environment;
    using namespace core::environment;

    //Tests that
    TorsionMoverOP allowed_mover( new TorsionMover( true, true ) );
    TorsionMoverOP duplicate_claim_mover( new TorsionMover( true, true ) );
    TorsionMoverOP no_claim_mover( new TorsionMover( false, true ) );
    TorsionMoverOP unreg_mover( new TorsionMover( true, true ) );

    EnvironmentOP env_op( new Environment( "torsion" ) );
    Environment & env = *env_op;

    env.register_mover( allowed_mover );
    env.register_mover( duplicate_claim_mover );
    env.register_mover( no_claim_mover );
    //don't register unreg_mover

    core::pose::Pose final_pose;

    {
      core::pose::Pose protected_pose = env.start( pose );
      // Verify conformation got copied into protected_pose pose.
      TS_ASSERT( protected_pose.total_residue() == pose.total_residue() );

      // Verify no_claim_mover can't change anything -- it shouldn't have a passport for this environment (NullPointer excn)
      TS_ASSERT_THROWS( no_claim_mover->apply( protected_pose ), EXCN_Env_Security_Exception );
      TS_ASSERT( protected_pose.phi( CLAIMED_RESID ) == init_phis[ CLAIMED_RESID ] );

      // Verify no_lock_mover can't change anything -- protected_pose shouldn't have a passport on its unlock stack
      TS_ASSERT_THROWS( allowed_mover->missing_unlock_apply( protected_pose ), EXCN_Env_Security_Exception );
      TS_ASSERT( protected_pose.phi( CLAIMED_RESID ) == init_phis[ CLAIMED_RESID ] );

      // Verify that unregistered mover lacks a passport for protected_pose conformation
      TS_ASSERT_THROWS( unreg_mover->apply( protected_pose ), utility::excn::EXCN_NullPointer );
      TS_ASSERT( protected_pose.phi( CLAIMED_RESID ) == init_phis[ CLAIMED_RESID ] );

      // Verify that allowed_mover can change it's claimed angle
      TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
      TS_ASSERT( protected_pose.phi( CLAIMED_RESID ) == NEW_PHI );

      // Verify that allowed_mover can do it again.
      TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
      TS_ASSERT( protected_pose.phi( CLAIMED_RESID ) == NEW_PHI );

      //Verify that duplicate mover *also* can make the change without throwing an error
      TS_ASSERT_THROWS_NOTHING( duplicate_claim_mover->apply( protected_pose ) );
      //Verify that allowed can't change 9 while claiming CLAIMED_RESID
      TS_ASSERT_THROWS( allowed_mover->apply( protected_pose, UNCLAIMED_RESID ), EXCN_Env_Security_Exception );
      TS_ASSERT_THROWS( duplicate_claim_mover->apply( protected_pose, UNCLAIMED_RESID ), EXCN_Env_Security_Exception );
      TS_ASSERT( protected_pose.phi( UNCLAIMED_RESID ) == init_phis[ UNCLAIMED_RESID ] );

      //Verify angles 1-9 are untouched in protected_pose
      for( core::Size i = 1; i <= pose.total_residue()-1; ++i ){
  TS_ASSERT_DELTA( pose.phi( i ), init_phis[i], 0.000001 );
      }

      TS_ASSERT_THROWS_NOTHING( final_pose = env.end( protected_pose ) );
    }

    //Verify angles 1-9 are untouched in pose and final_pose
    for( core::Size i = 1; i <= pose.total_residue(); ++i ){
      if( i != CLAIMED_RESID ){
  TS_ASSERT_DELTA( pose.phi( i ), final_pose.phi( i ), 0.000001 );
  TS_ASSERT_DELTA( pose.phi( i ), init_phis[i], 0.000001 );
      }
    }

    //Finish verification inital pose is unaffected
    TS_ASSERT_DELTA( pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 0.000001 );

    //Verify angle changes propagate out of environment
    TS_ASSERT_DELTA( final_pose.phi( CLAIMED_RESID ), NEW_PHI, 0.001 );

    //Verify other angle changes don't back-propagate to original pose
    TS_ASSERT_DIFFERS( pose.phi( CLAIMED_RESID ), final_pose.phi( CLAIMED_RESID ) );

    TS_TRACE( "End: test_single_phi_moves" );
  }

  void test_torsion_must_can_coexist(){
    TS_TRACE( "Beginning: test_torsion_must_can_coexist" );

    using namespace protocols::environment;
    using namespace core::environment;

    // Veirfy coexeistance of must and can control
    TorsionMoverOP must_mover( new TorsionMover( true, true, claims::MUST_CONTROL ) );
    TorsionMoverOP can_mover( new TorsionMover( true, true, claims::CAN_CONTROL ) );

    EnvironmentOP env_op( new Environment( "env" ) );
    Environment & env = *env_op;

    env.register_mover( must_mover );
    env.register_mover( can_mover );

    core::pose::Pose protected_pose;
    TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
    TS_ASSERT_THROWS_NOTHING( must_mover->apply( protected_pose ) );
    TS_ASSERT_THROWS_NOTHING( can_mover->apply( protected_pose ) );
    TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );

    TS_TRACE( "End: test_torsion_must_can_coexist" );
  }

  void test_torsion_can_exclusive_compatibility(){
    TS_TRACE( "Beginning: torsion_can_exclusive_compatibility" );

    using namespace protocols::environment;
    using namespace core::environment;

    TorsionMoverOP exclusive_mover( new TorsionMover( true, true, claims::EXCLUSIVE ) );
    TorsionMoverOP can_mover( new TorsionMover( true, true, claims::CAN_CONTROL ) );

    EnvironmentOP env_op( new Environment( "env" ) );
    Environment & env = *env_op;

    env.register_mover( exclusive_mover );
    env.register_mover( can_mover );

    core::pose::Pose protected_pose;
    TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
    TS_ASSERT_THROWS( can_mover->apply( protected_pose ), EXCN_Env_Security_Exception );
    TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ] );

    TS_ASSERT_THROWS_NOTHING( exclusive_mover->apply( protected_pose ) );
    TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );

    TS_TRACE( "End: torsion_can_exclusive_compatibility" );
  }

  void test_torsion_must_exclusive_incompatibility(){
    TS_TRACE( "Beginning: test_torsion_must_exclusive_incompatibility" );

    using namespace protocols::environment;
    using namespace core::environment;

    TorsionMoverOP exclusive_mover( new TorsionMover( true, true, claims::EXCLUSIVE ) );
    TorsionMoverOP must_mover( new TorsionMover( true, true, claims::MUST_CONTROL ) );

    EnvironmentOP env_op( new Environment( "env" ) );
    Environment & env = *env_op;

    env.register_mover( exclusive_mover );
    env.register_mover( must_mover );

    core::pose::Pose protected_pose;
    TS_ASSERT_THROWS( protected_pose = env.start( pose ), utility::excn::EXCN_BadInput );

    TS_TRACE( "End: test_torsion_must_exclusive_incompatibility" );
  }

  void test_torsion_init(){
    TS_TRACE( "Beginning: test_torsion_init" );

    using namespace protocols::environment;
    using namespace core::environment;

    TorsionMoverOP must_init_mover( new TorsionMover( true, true,
																											 claims::DOES_NOT_CONTROL,
																											 claims::EXCLUSIVE ) );
    TorsionMoverOP must_init_mover2( new TorsionMover( true, true,
																												claims::DOES_NOT_CONTROL,
																												claims::EXCLUSIVE ) );
    TorsionMoverOP active_can_init_mover( new TorsionMover( true, true,
																														 claims::DOES_NOT_CONTROL,
																														 claims::CAN_CONTROL ) );
    TorsionMoverOP inactive_can_init_mover( new TorsionMover( true, false,
																															 claims::DOES_NOT_CONTROL,
																															 claims::CAN_CONTROL ) );

    core::pose::Pose protected_pose;

    { // 2xMUST_INIT are not compatible
	  EnvironmentOP env_op( new Environment( "env" ) );
	  Environment & env = *env_op;

      env.register_mover( must_init_mover );
      env.register_mover( must_init_mover2 );
      TS_ASSERT_THROWS( protected_pose = env.start( pose ), utility::excn::EXCN_BadInput );
    }

    { // MUST_INIT takes precedence over CAN_INIT and does not complain.
      EnvironmentOP env_op( new Environment( "env" ) );
      Environment & env = *env_op;

      env.register_mover( must_init_mover );
      env.register_mover( inactive_can_init_mover );
      TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
      TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );
    }

    TS_TRACE( "End: test_torsion_init" );
  }

};
