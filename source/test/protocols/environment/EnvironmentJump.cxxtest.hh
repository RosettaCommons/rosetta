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
#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvClaimBroker.hh>

#include <test/protocols/environment/TestClaimingMover.hh>

//Other headers
#include <core/conformation/Conformation.hh>

#include <core/pack/task/residue_selector/ChainSelector.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/types.hh>

#include <test/core/init_util.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

#include <utility/tag/Tag.hh>

//C++ headers
#include <iostream>
#include <boost/bind/bind.hpp>

// --------------- Test Class --------------- //

const core::Size JUMP_START = 3;
const core::Size JUMP_END = 7;
const core::Size CUT_POS = 5;

class EnvironmentJump : public CxxTest::TestSuite {
public:

  // Shared data elements go here.
  core::pose::Pose pose;
  void (core::pose::Pose::*set_jump_) ( int const, const core::kinematics::Jump & );
	protocols::environment::claims::JumpClaimOP standard_claim_;

  // --------------- Fixtures --------------- //

  // Shared initialization goes here.
  void setUp() {
    core_init();

    using namespace protocols::environment;
    using namespace protocols::environment::claims;
    using namespace core::environment;

    core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

    set_jump_ = &core::pose::Pose::set_jump;
    standard_claim_ = new claims::JumpClaim( NULL, "claimed_jump",
                                             LocalPosition( "BASE", JUMP_START ),
                                             LocalPosition( "BASE", JUMP_END ) );
    standard_claim_->strength( claims::MUST_CONTROL, claims::DOES_NOT_CONTROL );
  }

  // Shared finalization goes here.
  void tearDown() {
  }

  void test_jump_moves() {
    TS_TRACE( "Beginning: test_jump_moves" );

    using namespace protocols::environment;
    using namespace protocols::environment::claims;
    using namespace core::environment;

    TesterOP allowed_mover = new Tester( new JumpClaim( *standard_claim_ ) );
    TesterOP duplicate_claim_mover = new Tester( new JumpClaim( *standard_claim_ ) );
    TesterOP no_claim_mover = new Tester();
    TesterOP unreg_mover = new Tester( new JumpClaim( *standard_claim_ ) );

    static_cast< JumpClaim* >( duplicate_claim_mover->claim().get() )->cut( core::environment::LocalPosition( "BASE", CUT_POS ) );

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

      core::kinematics::Jump new_jump( init_jump );
      new_jump.gaussian_move( 1, 10, 10 );

      //Verify invalid movers cannot move the jump, and do not change the conformation.
      TS_ASSERT_THROWS( no_claim_mover->apply( protected_pose, boost::bind( set_jump_, &protected_pose, 1, new_jump ) ),
                       EXCN_Env_Security_Exception );
      TS_ASSERT_THROWS( unreg_mover->apply( protected_pose, boost::bind( set_jump_, &protected_pose, 1, new_jump ) ),
                       utility::excn::EXCN_NullPointer );
      TS_ASSERT_EQUALS( protected_pose.jump( 1 ).get_translation(), init_jump.get_translation() );
      TS_ASSERT_EQUALS( protected_pose.jump( 1 ).get_rotation(), init_jump.get_rotation() );

      //Verify allowed movers are allowed
      TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose, boost::bind( set_jump_, &protected_pose, 1, new_jump ) ) );
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

  void test_autocut_placement( core::pose::Pose & pose ) {
    TS_TRACE( "Beginning: test_autocuts" );

    using namespace protocols::environment;
    using namespace protocols::environment::claims;
    using namespace core::environment;

    TesterOP allowed_mover = new Tester( new JumpClaim(*standard_claim_ ) );

    protocols::environment::Environment env( "env" );

    env.register_mover( allowed_mover );
    TS_ASSERT_THROWS( core::pose::Pose ppose = env.start( pose ), utility::excn::EXCN_BadInput );

    env.auto_cut( true );
    core::pose::Pose ppose;
    TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_jump(), 1 );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_cutpoint(), 1 );

    core::kinematics::Jump old_jump = ppose.jump( 1 );
    core::kinematics::Jump new_jump = old_jump;
    new_jump.gaussian_move( 1, 10.0, 10.0 );

    TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( ppose, boost::bind( set_jump_, &ppose, 1, new_jump ) ) );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ).get_translation(), old_jump.get_translation() );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ).get_rotation(), old_jump.get_rotation() );
    TS_ASSERT_THROWS_NOTHING( pose = env.end( ppose ) );

    TS_ASSERT_EQUALS( pose.fold_tree().num_jump(), 0 );
    TS_ASSERT_EQUALS( pose.fold_tree().num_cutpoint(), 0 );

    TS_TRACE( "End: test_autocuts" );
  }

  void test_cut_inheritance() {
    TS_TRACE( "Beginning: test_cut_inheritance" );

    using namespace protocols::environment;
		using namespace protocols::environment::claims;
    using namespace core::environment;

    core::environment::FoldTreeSketch new_fts( pose.total_residue() );
    new_fts.insert_cut( CUT_POS );
    new_fts.insert_jump( 1, pose.total_residue() );
    pose.fold_tree( *new_fts.render() );

    TesterOP allowed_mover = new Tester( new JumpClaim( *standard_claim_ ) );

    protocols::environment::Environment env( "env" );
    env.inherit_cuts( false );

    env.register_mover( allowed_mover );

    core::pose::Pose ppose;

    //Verify unbrokerable
    TS_ASSERT_THROWS( ppose = env.start( pose ) , utility::excn::EXCN_BadInput );

    //Verify allowing cut inheritance allows brokering.
    env.inherit_cuts( true );
    TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_jump(), 1 );
    TS_ASSERT_EQUALS( ppose.fold_tree().num_cutpoint(), 1 );
    TS_ASSERT( ppose.fold_tree().is_cutpoint( CUT_POS ) );

    //Verify jump exists, is movable.
    core::kinematics::Jump orig_jump = ppose.jump( 1 );
    core::kinematics::Jump new_jump( orig_jump );
    new_jump.gaussian_move( 1, 10.0, 10.0 );

    TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( ppose, boost::bind( set_jump_, &ppose, 1, new_jump ) ) );
    TS_ASSERT_DIFFERS( ppose.jump( 1 ), orig_jump );

    //Verify pose closable.
    core::pose::Pose final_pose;
    TS_ASSERT_THROWS_NOTHING( final_pose = env.end( ppose ) );

    TS_ASSERT_EQUALS( final_pose.fold_tree(), pose.fold_tree() );

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
    core::pose::PDBInfoOP info = new core::pose::PDBInfo( pose.total_residue() );
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

  void test_JumpClaim_rscripts_creation( core::pose::Pose const& pose ) {
    using namespace core::pack::task::residue_selector;
    using namespace protocols::environment;
    using namespace protocols::environment::claims;

    core::pack::task::residue_selector::ChainSelectorOP chA_sele = new core::pack::task::residue_selector::ChainSelector();
    chA_sele->set_chain_strings( utility::vector1< std::string >( 1, "1" ) );

    core::pack::task::residue_selector::ChainSelectorOP chB_sele = new core::pack::task::residue_selector::ChainSelector();
    chB_sele->set_chain_strings( utility::vector1< std::string >( 1, "2" ) );

    basic::datacache::DataMap datamap;
    datamap.add( "ResidueSelector", "ChainA", chA_sele );
    datamap.add( "ResidueSelector", "ChainB", chB_sele );

    std::string const tag_string = "<JumpClaim jump_label=\"labelA\" control_strength=CAN_CONTROL position1=\"ChainA,5\" position2=\"ChainB,5\" />";
    std::stringstream ss( tag_string );
    utility::tag::TagPtr tag = new utility::tag::Tag;
    tag->read( ss );

    TesterOP tester = new Tester();
    tester->claim( EnvClaim::make_claim( tag->getName(), tester, tag, datamap ) );
    protocols::environment::Environment env( "test" );
    env.register_mover( tester );

    core::pose::Pose ppose;
    TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );

    core::kinematics::Jump const old_jump = ppose.jump( 1 );
    core::kinematics::Jump new_jump( old_jump );
    new_jump.gaussian_move( 1, 2.0, 2.0 );

    TS_ASSERT_THROWS_NOTHING( tester->apply( ppose, boost::bind( boost::bind( set_jump_, &ppose, 1, new_jump ) ) ) );

  }
};
