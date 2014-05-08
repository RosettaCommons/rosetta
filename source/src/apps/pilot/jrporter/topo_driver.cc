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

#include <core/id/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/environment/DofPassport.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/abinitio/abscript/RigidChunkCM.hh>

#include <protocols/simple_moves/FragmentMover.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

static basic::Tracer tr("main");

int main( int argc, char** argv ){

  using namespace protocols;
  using namespace protocols::environment;


  devel::init( argc, argv );

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, "FRIENDLYFRIENDSSSSSSSSSSSS", "fa_standard");
  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    pose.set_phi( i, 180 );
    pose.set_psi( i, 180 );
    pose.set_omega( i, 180 );
  }

  core::pose::Pose templ;
  core::pose::make_pose_from_sequence(templ, "FRIENDLYFRIENDSSSSSSSSSSSS", "fa_standard");
  for( core::Size i = 1; i <= templ.total_residue(); ++i ){
    templ.set_phi( i, -64.0 );
    templ.set_psi( i, -41.0 );
    templ.set_omega( i, 180 );
  }

  std::pair< core::Size, core::Size > l1 = std::make_pair( 3, 6 );
  std::pair< core::Size, core::Size > l2 = std::make_pair( 12, 15 );

  loops::SerializedLoopList l;
  loops::SerializedLoop p;
  p.start = l1.first; p.stop = l1.second; p.skip_rate = 0.0; p.extended = false;
  l.push_back( p );
  p.start = l2.first; p.stop = l2.second; p.skip_rate = 0.0; p.extended = false;
  l.push_back( p );


  abinitio::abscript::RigidChunkCMOP m = new abinitio::abscript::RigidChunkCM( "BASE",
                                                                               loops::Loops( l ),
                                                                               templ );

  core::pose::Pose end;
  core::pose::Pose ppose;
  {
    Environment env( "env" );
    env.register_mover( ClaimingMoverOP( m ) );

    ppose = env.start( pose );

    std::cout << ppose.fold_tree() << std::endl;
    assert( ppose.fold_tree().num_jump() == 1 );
    assert( ppose.fold_tree().cutpoint(1) > (int) l1.second &&
            ppose.fold_tree().cutpoint(1) < (int) l2.first );

    end = env.end( ppose );
  }

  templ.dump_pdb( "tmpl.pdb" );
  end.dump_pdb( "tmp.pdb" );

  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    using namespace core::id;

    if( ( i >= l1.first && i <= l1.second ) ||
        ( i >= l2.first && i <= l2.second ) ){
      // inside chunk
      if( !ppose.fold_tree().is_cutpoint( (int) i-1 ) && i > 0 ){
        assert( std::abs( end.phi(i) - templ.phi(i) ) < 0.000001 );
      }
      if( !ppose.fold_tree().is_cutpoint( (int) i ) ){
      assert( std::abs( end.psi(i) - templ.psi(i) ) < 0.000001 );
      }
    } else  if( i == l1.first - 1 ||
                i == l2.first - 1 ) {
      // residues before the chunk are a mix
      if( !ppose.fold_tree().is_cutpoint( (int) i-1 ) && i > 0 ){
        assert( std::abs( end.phi(i) - pose.phi(i) ) < 0.000001 );
      }
      if( !ppose.fold_tree().is_cutpoint( (int) i ) ){
        assert( std::abs( end.psi(i) - templ.psi(i) ) < 0.000001 );
      }

      // require residue connections to the chunk
      assert( end.conformation().dof_id_from_torsion_id( TorsionID( i, BB, psi_torsion ) ).valid() );
    } else if ( i == l1.second + 1 ||
                i == l2.second + 1) {
      // residues after the chunk are a mix
      if( end.conformation().dof_id_from_torsion_id( TorsionID( i, BB, phi_torsion ) ).valid() ){
        assert( std::abs( end.phi(i) - templ.phi(i) ) < 0.000001 );
      }

      if( i < ppose.total_residue() && !ppose.fold_tree().is_cutpoint( (int) i ) ){
        assert( std::cos( end.psi(i) ) - std::cos( pose.psi(i) ) < 0.000001 );
      }

      //require residue connections to the chunk
      assert( end.conformation().dof_id_from_torsion_id( TorsionID( i, BB, phi_torsion ) ).valid() );
    } else {
      // outside chunk
      if( end.conformation().dof_id_from_torsion_id( TorsionID( i, BB, phi_torsion ) ).valid() ){
        assert( ( end.phi(i) - pose.phi(i) ) < 0.000001 );
      }
      if( end.conformation().dof_id_from_torsion_id( TorsionID( i, BB, psi_torsion ) ).valid() ){
        assert( ( end.psi(i) - pose.psi(i) ) < 0.000001 );
      }
    }
  }

  tr.Info << "Done!" << std::endl;
}