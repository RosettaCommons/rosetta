// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file topo_driver.cc
/// @brief  Sandbox executable for the Environment and associated classes.
/// @author Justin Porter

#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/abinitio/abscript/RigidChunkCM.hh>
#include <protocols/abinitio/abscript/FragmentCM.hh>

#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/EnvMover.hh>
#include <protocols/environment/Environment.hh>

#include <test/protocols/environment/TorsionMover.hh>

#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>

#include <devel/init.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

#include <core/id/types.hh>

#include <boost/bind/bind.hpp>

#define TS_ASSERT( x ) assert( x )
#define TS_ASSERT_THROWS_NOTHING( x ) x
#define TS_ASSERT_EQUALS( x , y ) TS_ASSERT( x == y )
#define TS_ASSERT_DIFFERS( x , y ) TS_ASSERT( x != y )
#define TS_ASSERT_DELTA( x , y, d ) TS_ASSERT( std::abs( x - y ) < d )
#define TS_ASSERT_THROWS( x , y ) try{ x; } catch( y ){}
#define TS_ASSERT_LESS_THAN( x, y ) TS_ASSERT( x < y );
#define TS_TRACE( x ) std::cout << x << std::endl;
#define ANGLE_DELTA( x , y , d ) TS_ASSERT_DELTA( std::cos( x ) , std::cos( y ), d );

std::string const FRAGFILE_LOCATION = "/Users/jrporter/Rosetta/main/source/test/protocols/abinitio/abscript/one_frag3_per_pos";

void test_chunk( core::pose::Pose const pose );

void test_chunk( core::pose::Pose const pose ) {
  
  utility::vector1< core::Real > init_phis;
  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    init_phis.push_back( pose.phi( i ) );
  }


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
    TS_ASSERT_EQUALS( protected_pose.total_residue(), pose.total_residue() );

    // Verify no_claim_mover can't change anything -- it shouldn't have a passport for this environment (NullPointer excn)
    TS_ASSERT_THROWS( no_claim_mover->apply( protected_pose ), EXCN_Env_Security_Exception );
    TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

    // Verify no_lock_mover can't change anything -- protected_pose shouldn't have a passport on its unlock stack
    TS_ASSERT_THROWS( allowed_mover->missing_unlock_apply( protected_pose ), EXCN_Env_Security_Exception );
    TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

    // Verify that unregistered mover lacks a passport for protected_pose conformation
    TS_ASSERT_THROWS( unreg_mover->apply( protected_pose ), utility::excn::EXCN_NullPointer );
    TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

    // Verify that allowed_mover can change it's claimed angle
    TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
    TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), NEW_PHI, 1e-12 );

    // Verify that allowed_mover can do it again.
    TS_ASSERT_THROWS_NOTHING( allowed_mover->apply( protected_pose ) );
    TS_ASSERT_DELTA( protected_pose.phi( CLAIMED_RESID ), NEW_PHI, 1e-10 );

    //Verify that duplicate mover *also* can make the change without throwing an error
    TS_ASSERT_THROWS_NOTHING( duplicate_claim_mover->apply( protected_pose ) );
    //Verify that allowed can't change 9 while claiming CLAIMED_RESID
    TS_ASSERT_THROWS( allowed_mover->apply( protected_pose, UNCLAIMED_RESID ), EXCN_Env_Security_Exception );
    TS_ASSERT_THROWS( duplicate_claim_mover->apply( protected_pose, UNCLAIMED_RESID ), EXCN_Env_Security_Exception );
    TS_ASSERT_DELTA( protected_pose.phi( UNCLAIMED_RESID ), init_phis[ UNCLAIMED_RESID ], 1e-12 );

    //Verify angles 1-9 are untouched in protected_pose
    for( core::Size i = 1; i <= pose.total_residue()-1; ++i ){
      TS_ASSERT_DELTA( pose.phi( i ), init_phis[i], 1e-12 );
    }

    TS_ASSERT_THROWS_NOTHING( final_pose = env.end( protected_pose ) );
  }

  // Verify angles 1-9 are untouched in pose and final_pose;
  // Phi 1 is not well-defined, so skip that one.
  for( core::Size i = 2; i <= pose.total_residue(); ++i ){
    if( i != CLAIMED_RESID ){
      TS_ASSERT_DELTA( pose.phi( i ), final_pose.phi( i ), 1e-12 );
      TS_ASSERT_DELTA( pose.phi( i ), init_phis[i], 1e-12 );
    }
  }

  //Finish verification inital pose is unaffected
  TS_ASSERT_DELTA( pose.phi( CLAIMED_RESID ), init_phis[ CLAIMED_RESID ], 1e-12 );

  //Verify angle changes propagate out of environment
  TS_ASSERT_DELTA( final_pose.phi( CLAIMED_RESID ), NEW_PHI, 0.001 );

  //Verify other angle changes don't back-propagate to original pose
  TS_ASSERT_DIFFERS( pose.phi( CLAIMED_RESID ), final_pose.phi( CLAIMED_RESID ) );

  TS_TRACE( "End: test_single_phi_moves" );
}

int main( int argc, char** argv ){
  devel::init( argc, argv );

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, "FRMQIFVYFRIENDS", core::chemical::FA_STANDARD);

  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    pose.set_phi( i, -65 );
    pose.set_psi( i, -41 );
    pose.set_omega( i, 180 );
  }

  try {
    test_chunk( pose );
  } catch ( utility::excn::EXCN_Msg_Exception excn ){
    std::cout << excn << std::endl;
    std::exit( 1 );
  }

  std::cout << "Success!" << std::endl;

}
