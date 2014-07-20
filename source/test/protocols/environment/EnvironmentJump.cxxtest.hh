// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

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
#include <core/pose/PDB_Info.hh>
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

const core::Size JUMP_START = 3;
const core::Size JUMP_END = 7;
const core::Size CUT_POS = 5;

class ToyMover : public protocols::environment::ClaimingMover {
protected:

  ToyMover( bool move ):
  move_( move )
  {}

  //Yes I know this is against the law
  bool move_;
};

class JumpMover : public ToyMover {
public:
  claims::JumpClaimOP claim_;

  JumpMover( bool claim, bool move ):
  ToyMover( move ),
  claim_( 0 )
  {

    using core::environment::LocalPosition;

    if( claim ){
      claim_ = new claims::JumpClaim( this,
                                     "claimed_jump",
                                     LocalPosition( "BASE", JUMP_START ),
                                     LocalPosition( "BASE", JUMP_END ) );
      claim_->strength( claims::MUST_CONTROL, claims::DOES_NOT_CONTROL );
    }

  }

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

  claims::JumpClaimOP claim() {
    return claim_;
  }

  virtual claims::EnvClaims yield_claims( core::pose::Pose const&,
                                         basic::datacache::WriteableCacheableMapOP ){
    claims::EnvClaims claims;
    if( claim_ ){
      claims.push_back( claim_ );
    }
    return claims;
  }
};

typedef utility::pointer::owning_ptr< JumpMover > JumpMoverOP;

} //environment
} //protocols

// --------------- Test Class --------------- //

class EnvironmentJump : public CxxTest::TestSuite {
public:

  // Shared data elements go here.
  core::pose::Pose pose;

  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
    core_init();

    using namespace protocols::environment;
    using namespace core::environment;

    core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_jump_moves() {
    TS_TRACE( "Beginning: test_jump_moves" );

    using namespace protocols::environment;
    using namespace core::environment;

    JumpMoverOP allowed_mover = new JumpMover( true, true );
    JumpMoverOP duplicate_claim_mover = new JumpMover( true, true );
    JumpMoverOP no_claim_mover = new JumpMover( false, true );
    JumpMoverOP unreg_mover = new JumpMover( true, true );

    duplicate_claim_mover->claim()->cut( core::environment::LocalPosition( "BASE", CUT_POS ) );

    protocols::environment::Environment env( "env" );

    env.auto_cut( true );

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
      TS_ASSERT( protected_pose.fold_tree().is_cutpoint( CUT_POS ) );

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

      TS_ASSERT_LESS_THAN( ( (protected_pose.jump( 1 ).get_translation() - re_ft_app.jump( 1 ).get_translation() ).length() ), 0.000001 );

      core::Real rotation_delta = ( ( protected_pose.jump( 1 ).get_rotation() - re_ft_app.jump( 1 ).get_rotation() ).row_x().length() ) +
                                  ( ( protected_pose.jump( 1 ).get_rotation() - re_ft_app.jump( 1 ).get_rotation() ).row_y().length() ) +
                                  ( ( protected_pose.jump( 1 ).get_rotation() - re_ft_app.jump( 1 ).get_rotation() ).row_z().length() );

      TS_ASSERT_LESS_THAN( rotation_delta, 0.000001 );
    }

    TS_TRACE( "End: test_jump_moves" );
  }

  void test_autocut_placement() {
    TS_TRACE( "Beginning: test_autocuts" );

    using namespace protocols::environment;
    using namespace core::environment;

    JumpMoverOP allowed_mover = new JumpMover( true, true );

    protocols::environment::Environment env( "env" );

    env.register_mover( allowed_mover );
    TS_ASSERT_THROWS( core::pose::Pose ppose = env.start( pose ), utility::excn::EXCN_BadInput );

    env.auto_cut( true );
    core::pose::Pose ppose;
    TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_jump(), 1 );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_cutpoint(), 1 );
    core::kinematics::Jump orig_jump = ppose.jump( 1 );
    TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( ppose ) );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ).get_translation(), orig_jump.get_translation() );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ).get_rotation(), orig_jump.get_rotation() );
    TS_ASSERT_THROWS_NOTHING( pose = env.end( ppose ) );

    TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 0 );
    TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 0 );

    TS_TRACE( "End: test_autocuts" );
  }

  void test_cut_inheritance() {
    TS_TRACE( "Beginning: test_cut_inheritance" );

    using namespace protocols::environment;
    using namespace core::environment;

    core::environment::FoldTreeSketch new_fts( pose.total_residue() );
    new_fts.insert_cut( CUT_POS );
    new_fts.insert_jump( 1, pose.total_residue() );
    core::kinematics::FoldTree new_ft( *new_fts.render() );
    pose.fold_tree( new_ft );

    JumpMoverOP allowed_mover = new JumpMover( true, true );

    protocols::environment::Environment env( "env" );
    env.inherit_cuts( false );

    env.register_mover( allowed_mover );

    core::pose::Pose ppose;
    TS_ASSERT_THROWS( ppose = env.start( pose ) , utility::excn::EXCN_BadInput );

    env.inherit_cuts( true );
    TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_jump(), 1 );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_cutpoint(), 1 );
    TS_ASSERT( ppose.fold_tree().is_cutpoint( CUT_POS ) );
    core::kinematics::Jump orig_jump = ppose.jump( 1 );
    TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( ppose ) );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ).get_translation(), orig_jump.get_translation() );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ).get_rotation(), orig_jump.get_rotation() );
    TS_ASSERT_THROWS_NOTHING( pose = env.end( ppose ) );

    TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 1 );
    TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 1 );
    TS_ASSERT_EQUALS( pose.fold_tree().upstream_jump_residue( 1 ), 1 );
    TS_ASSERT_EQUALS( pose.fold_tree().downstream_jump_residue( 1 ), pose.total_residue() );

    TS_TRACE( "End: test_cut_inheritance" );
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

  void test_pdb_info_persistence() {
    core::pose::PDB_InfoOP info = new core::pose::PDB_Info( pose.total_residue() );
    info->set_chains( 'A' );
    pose.pdb_info( info );

    protocols::environment::Environment env( "env" );

    core::pose::Pose protected_pose;
    TS_ASSERT_THROWS_NOTHING( protected_pose = env.start( pose ) );
    TS_ASSERT( protected_pose.pdb_info() );
    TS_ASSERT_EQUALS( protected_pose.pdb_info()->chain( 1 ), 'A' );

    core::pose::Pose end_pose;
    TS_ASSERT_THROWS_NOTHING( end_pose = env.end( protected_pose ) );
    TS_ASSERT( end_pose.pdb_info() );
    TS_ASSERT_EQUALS( end_pose.pdb_info()->chain( 1 ), 'A' );
  }
};
