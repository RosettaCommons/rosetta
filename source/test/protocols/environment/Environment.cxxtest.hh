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
#include <protocols/environment/ClaimingMover.hh>
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

// ---------------- Toy Movers --------------- //

namespace protocols {
namespace environment {

const core::Size CLAIMED_RESID = 10;
const core::Size UNCLAIMED_RESID = 9;
const core::Real NEW_PHI = 23.0;

const core::Size JUMP_START = 3;
const core::Size JUMP_END = 7;
class ToyMover : public protocols::environment::ClaimingMover {
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
    claims::InitializationStrength i_str = claims::DOES_NOT_INITIALIZE):
    ToyMover( claim, move ),
    resnum_( CLAIMED_RESID ),
    control_str_( c_str ),
    init_str_( i_str )
  {}

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
      claims::TorsionClaimOP new_claim = new claims::TorsionClaim( this, LocalPosition( "BASE", resnum_ ) );
      new_claim->ctrl_strength( control_str_ );
      new_claim->init_strength( init_str_ );
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
  claims::InitializationStrength init_str_;
};

class JumpMover : public ToyMover {
public:
  JumpMover( bool claim, bool move ):
    ToyMover( claim, move )
  {}

  virtual void apply( Pose& pose ){
    apply( pose, 1 );
  }

  virtual void apply( Pose& pose, Size jumpnum ){
    DofUnlock activation( pose.conformation(), passport() );
    core::kinematics::Jump j;
    j.set_translation( numeric::xyzVector< core::Real >( 5, 10, 15 ) );
    if( move_ ){
      pose.set_jump( jumpnum, j );
    }
  }

  virtual std::string get_name() const{
    return "JumpMover";
  }

  virtual claims::EnvClaims yield_claims( core::pose::Pose const&,
					  basic::datacache::WriteableCacheableMapOP ){
    using core::environment::LocalPosition;
    claims::EnvClaims claims;

    if( claim_ ){
      claims::JumpClaimOP c = new claims::JumpClaim( this,
                 "claimed_jump",
                 LocalPosition( "BASE", 3 ),
                 LocalPosition( "BASE", 7 ) );
      c->ctrl_strength( claims::MUST_CONTROL );
      claims.push_back( c );
    }

    return claims;
  }
};


typedef utility::pointer::owning_ptr< TorsionMover > TorsionMoverOP;
typedef utility::pointer::owning_ptr< JumpMover > JumpMoverOP;

} //environment
} //protocols

// --------------- Test Class --------------- //

class Environment : public CxxTest::TestSuite {
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

    TorsionMoverOP mover = new TorsionMover( false, false );

    protocols::environment::Environment env( "env" );

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

    TorsionMoverOP tier1_mover = new TorsionMover( true, true );
    TorsionMoverOP tier2_mover = new TorsionMover( true, true );

    protocols::environment::Environment env1( "tier1" );
    protocols::environment::Environment env2( "tier2" );

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
    TorsionMoverOP allowed_mover = new TorsionMover( true, true );
    TorsionMoverOP duplicate_claim_mover = new TorsionMover( true, true );
    TorsionMoverOP no_claim_mover = new TorsionMover( false, true );
    TorsionMoverOP unreg_mover = new TorsionMover( true, true );

    protocols::environment::Environment env( "torsion" );

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
    TorsionMoverOP must_mover = new TorsionMover( true, true, claims::MUST_CONTROL );
    TorsionMoverOP can_mover = new TorsionMover( true, true, claims::CAN_CONTROL );

    protocols::environment::Environment env( "env" );

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

    TorsionMoverOP exclusive_mover = new TorsionMover( true, true, claims::EXCLUSIVE );
    TorsionMoverOP can_mover = new TorsionMover( true, true, claims::CAN_CONTROL );

    protocols::environment::Environment env( "env" );

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

    TorsionMoverOP exclusive_mover = new TorsionMover( true, true, claims::EXCLUSIVE );
    TorsionMoverOP must_mover = new TorsionMover( true, true, claims::MUST_CONTROL );

    protocols::environment::Environment env( "env" );

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

    TorsionMoverOP must_init_mover = new TorsionMover( true, true,
                   claims::DOES_NOT_CONTROL,
                   claims::MUST_INITIALIZE );
    TorsionMoverOP must_init_mover2 = new TorsionMover( true, true,
              claims::DOES_NOT_CONTROL,
              claims::MUST_INITIALIZE );
    TorsionMoverOP active_can_init_mover = new TorsionMover( true, true,
                   claims::DOES_NOT_CONTROL,
                   claims::CAN_INITIALIZE);
    TorsionMoverOP inactive_can_init_mover = new TorsionMover( true, false,
                   claims::DOES_NOT_CONTROL,
                   claims::CAN_INITIALIZE);

    core::pose::Pose protected_pose;

    { // 2xMUST_INIT are not compatible
      protocols::environment::Environment env( "env" );

      env.register_mover( must_init_mover );
      env.register_mover( must_init_mover2 );
      TS_ASSERT_THROWS( protected_pose = env.start( pose ), utility::excn::EXCN_BadInput );
    }

    { // MUST_INIT takes precedence over CAN_INIT and does not complain.
      protocols::environment::Environment env( "env" );

      env.register_mover( must_init_mover );
      env.register_mover( inactive_can_init_mover );
      TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
      TS_ASSERT_EQUALS( protected_pose.phi( CLAIMED_RESID ), NEW_PHI );
    }

    TS_TRACE( "End: test_torsion_init" );
  }

  void test_jump_moves() {
    TS_TRACE( "Beginning: test_jump_moves" );

    using namespace protocols::environment;
    using namespace core::environment;

    //Tests that
    JumpMoverOP allowed_mover = new JumpMover( true, true );
    JumpMoverOP duplicate_claim_mover = new JumpMover( true, true );
    JumpMoverOP no_claim_mover = new JumpMover( false, true );
    JumpMoverOP unreg_mover = new JumpMover( true, true );

    protocols::environment::Environment env( "env" );

    env.register_mover( allowed_mover );
    env.register_mover( duplicate_claim_mover );
    env.register_mover( no_claim_mover );
    // do not register unreg_mover

    core::pose::Pose final_pose;

    {
      core::pose::Pose protected_pose;
      TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
      core::kinematics::FoldTree prot_ft ( protected_pose.fold_tree() );

      TS_ASSERT_EQUALS( protected_pose.sequence(), pose.sequence() );
      TS_ASSERT_DIFFERS( protected_pose.annotated_sequence(), pose.annotated_sequence() );

      // Verify jump now exists
      core::kinematics::Jump init_jump;
      TS_ASSERT_EQUALS( protected_pose.num_jump(), core::Size( 1 ) );
      TS_ASSERT_THROWS_NOTHING( init_jump = protected_pose.jump( 1 ) );

      //Verify invalid movers cannot move the jump, and do not change the conformation.
      TS_ASSERT_THROWS( no_claim_mover->apply( protected_pose ), EXCN_Env_Security_Exception );
      TS_ASSERT_THROWS( unreg_mover->apply( protected_pose ), utility::excn::EXCN_NullPointer );
      TS_ASSERT_EQUALS( protected_pose.jump( 1 ).get_translation(), init_jump.get_translation() );
      TS_ASSERT_EQUALS( protected_pose.jump( 1 ).get_rotation(), init_jump.get_rotation() );

      //Verify allowed movers are allowed
      TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
      TS_ASSERT_DIFFERS( protected_pose.jump( 1 ).get_translation(), init_jump.get_translation() );
      TS_ASSERT_DIFFERS( protected_pose.jump( 1 ).get_rotation(), init_jump.get_rotation() );

      TS_ASSERT_THROWS_NOTHING( final_pose = env.end( protected_pose ) );

      // Verify fold tree reset.
      TS_ASSERT_EQUALS( final_pose.fold_tree().num_cutpoint(), 0 );

      //Verify jump RT changes are passed through
      core::pose::Pose re_ft_app( final_pose );
      re_ft_app.fold_tree( prot_ft );
      final_pose.dump_pdb( "finalpose.pdb" );
      protected_pose.dump_pdb( "protectedpose.pdb" );

      TS_ASSERT_LESS_THAN( ( (protected_pose.jump( 1 ).get_translation() - re_ft_app.jump( 1 ).get_translation() ).length() ), 0.000001 );

      core::Real rotation_delta = ( ( protected_pose.jump( 1 ).get_rotation() - re_ft_app.jump( 1 ).get_rotation() ).row_x().length() ) +
	( ( protected_pose.jump( 1 ).get_rotation() - re_ft_app.jump( 1 ).get_rotation() ).row_y().length() ) +
	( ( protected_pose.jump( 1 ).get_rotation() - re_ft_app.jump( 1 ).get_rotation() ).row_z().length() );

      TS_ASSERT_LESS_THAN( rotation_delta, 0.000001 );
    }

    TS_TRACE( "End: test_jump_moves" );
  }

  void test_cache_persistence() {
    TS_TRACE( "Beginning: test_cache_persistence" );

    using namespace core::pose::datacache;
    using namespace basic::datacache;

    pose.data().set( CacheableDataType::WRITEABLE_DATA, new WriteableCacheableMap() );
    core::pose::Pose protected_pose;
    core::pose::Pose final_pose;

    protocols::environment::Environment env( "env" );

    TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
    TS_ASSERT_THROWS_NOTHING( final_pose = env.end( protected_pose ) );

    final_pose.data().has( CacheableDataType::WRITEABLE_DATA );

    TS_ASSERT( final_pose.data().has( CacheableDataType::WRITEABLE_DATA ) );

    TS_TRACE( "End: test_jump_moves" );
  }
};
