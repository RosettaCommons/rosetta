// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file topo_driver.cc
/// @brief  Sandbox executable for the Environment and associated classes.
/// @author Justin Porter

#include <basic/Tracer.hh>


#include <basic/Tracer.hh>

#include <core/id/types.hh>

#include <core/pack/task/residue_selector/ChainSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/kinematics/Jump.hh>

#include <protocols/environment/Environment.hh>
#include <protocols/environment/CoMTrackerCM.hh>
#include <protocols/environment/EnvExcn.hh>

#include <protocols/environment/claims/JumpClaim.hh>

#include <protocols/abinitio/abscript/StructPerturberCM.hh>

#include <protocols/rigid/UniformRigidBodyCM.hh>

#include <devel/init.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#define TS_ASSERT( x ) assert( x )
#define TS_ASSERT_THROWS_NOTHING( x ) x
#define TS_ASSERT_EQUALS( x , y ) TS_ASSERT( x == y )
#define TS_ASSERT_DIFFERS( x , y ) TS_ASSERT( x != y )
#define TS_ASSERT_DELTA( x , y, d ) TS_ASSERT( std::abs( x - y ) < d )
#define TS_ASSERT_THROWS( x , y ) try{ x; } catch( y ){}
#define TS_ASSERT_LESS_THAN( x, y ) TS_ASSERT( x < y );
#define TS_TRACE( x ) std::cout << x << std::endl;

static basic::Tracer tr("main");

void test_com_rb( core::pose::Pose pose ) {

  using namespace core::environment;
  using namespace protocols;
  using namespace protocols::environment;
  using namespace core::pack::task::residue_selector;
  using namespace numeric;
  using namespace core::conformation;
  using namespace abinitio::abscript;
  using namespace protocols::rigid;

  TS_ASSERT_DIFFERS( pose.total_residue(), 0 );

  std::string const com1 = "com1";
  std::string const com2 = "com2";

  ChainSelectorOP chain1 = new ChainSelector();
  ChainSelectorOP chain2 = new ChainSelector();
  chain1->set_chain_strings( utility::vector1< std::string >( 1, "1" ) );
  chain2->set_chain_strings( utility::vector1< std::string >( 1, "2" ) );
  CoMTrackerCMOP tracker1 = new CoMTrackerCM( com1, chain1 );
  CoMTrackerCMOP tracker2 = new CoMTrackerCM( com2, chain2 );

  StructPerturberCMOP perturb = new StructPerturberCM( "BASE", 10.0 );

  UniformRigidBodyCMOP rigpert = new UniformRigidBodyCM( "perturb", LocalPosition( com1, 1 ), LocalPosition( com2, 1 ) );

  Environment env( "env" );
  env.register_mover( rigpert );
  env.register_mover( tracker1 );
  env.register_mover( tracker2 );
  env.register_mover( perturb );

  core::pose::Pose ppose;

  //Test correct topology construction
  TS_ASSERT_THROWS_NOTHING( ppose = env.start( pose ) );
  TS_ASSERT_EQUALS( ppose.num_jump(), 3 );

  // TODO: make this determination less hacky and more declarative
  int const comcom_jump = (int) ppose.fold_tree().jump_nr( ppose.total_residue()-1, ppose.total_residue() );
  int const chain1_com =  (int) ppose.fold_tree().jump_nr( 2, ppose.total_residue()-1 );
  int const chain2_com =  (int) ppose.fold_tree().jump_nr( ppose.total_residue()/2+1, ppose.total_residue() );

  TS_ASSERT_DIFFERS( comcom_jump, 0 );
  TS_ASSERT_DIFFERS( chain1_com, 0 );
  TS_ASSERT_DIFFERS( chain2_com, 0 );

  //Verify initialization perturbations
  utility::vector1< core::kinematics::Jump > jumps( 3 );
  for( Size i = 1; i <= ppose.num_jump(); ++i ){
    core::kinematics::Jump j;
    {
      // calculate original transform between com1 and
      core::pose::Pose tmp( pose );
      core::kinematics::FoldTree tmp_ft = ppose.fold_tree();
      tmp_ft.delete_seqpos( (int) pose.total_residue() );
      tmp_ft.delete_seqpos( (int) pose.total_residue()-1 );
      tmp.fold_tree( tmp_ft );
      j = tmp.jump( (int) i );
    }

    TS_ASSERT_DIFFERS( j.get_translation(), ppose.jump( (int) i ).get_translation() );
    TS_ASSERT_DIFFERS( j.get_rotation(), ppose.jump( (int) i ).get_rotation() );

    //store values of jumps for later comparison
    jumps[i] = ppose.jump( (int) i );
  }

  //Verify rigpert perturbation.
  TS_ASSERT_THROWS_NOTHING( rigpert->apply( ppose ) );

  // Verify rigpert modifies no jumps
  TS_ASSERT_DIFFERS( jumps[comcom_jump].get_translation(), ppose.jump( comcom_jump ).get_translation() );
  TS_ASSERT_DIFFERS( jumps[comcom_jump].get_rotation(), ppose.jump( comcom_jump ).get_rotation() );
  TS_ASSERT_LESS_THAN( ( jumps[chain1_com].get_translation() - ppose.jump( chain1_com ).get_translation() ).length(), 1e-10 );
  TS_ASSERT_LESS_THAN( ( jumps[chain2_com].get_translation() - ppose.jump( chain2_com ).get_translation() ).length(), 1e-10 );

  // Recalculate CoM residues.
  TS_ASSERT_THROWS_NOTHING( tracker1->apply( ppose ) );
  TS_ASSERT_THROWS_NOTHING( tracker2->apply( ppose ) );

  // Verify correct recalculation of CoM Residues
  TS_ASSERT_DIFFERS( jumps[comcom_jump].get_translation(), ppose.jump( comcom_jump ).get_translation() );
  TS_ASSERT_DIFFERS( jumps[comcom_jump].get_rotation(), ppose.jump( comcom_jump ).get_rotation() );
  TS_ASSERT_LESS_THAN( ( jumps[chain1_com].get_translation() - ppose.jump( chain1_com ).get_translation() ).length(), 1e-10 );
  TS_ASSERT_LESS_THAN( ( jumps[chain2_com].get_translation() - ppose.jump( chain2_com ).get_translation() ).length(), 1e-10 );

  //Test environment end
  core::pose::Pose final_pose;
  TS_ASSERT_THROWS_NOTHING( final_pose = env.end( ppose ) );
}

int main( int argc, char** argv ){
  devel::init( argc, argv );

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDS", "fa_standard");

  for( Size i = 1; i <= pose.total_residue(); ++i ){
    pose.set_phi( i, -65 );
    pose.set_psi( i, -41 );
  }

  pose.append_pose_by_jump( *pose.clone(), 5 );
  core::kinematics::Jump j = pose.jump( 1 );
  j.gaussian_move(1, 50, 0);
  pose.set_jump(1, j);

  try {
    test_com_rb( pose );
  } catch ( utility::excn::EXCN_Msg_Exception excn ){
    std::cout << excn << std::endl;
    std::exit( 1 );
  }

  std::cout << "Done!" << std::endl;

}
