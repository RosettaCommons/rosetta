// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-

//Auto Headers
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KinematicTaskCenter
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @detailed responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange

namespace protocols {


core::pose::Pose ResolutionSwitcher::start_pose() {
  if ( start_centroid_ ) {
    if ( init_fa_ ) {
      pose::Pose pose( init_pose_ );
      core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
      return pose;
    } else return init_pose_;
  } else if ( init_fa_ ) {
    return init_pose_;
  } else {
    utility_exit_with_message("don't have full-atom pose to start");
    return init_pose_; //to make compiler happy.
    //or switch to full-atom?
  }
}



bool copy_side_chains(
   core::pose::Pose& pose,
   utility::vector1< bool >& needToRepack,
   core::pose::Pose const& fa_input_pose
) {
  // copy sidechain torsions from input pose
  assert( pose.total_residue() == fa_input_pose.total_residue() );

  for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
    // if there is missing density in the sidechain, then we need to repack
    // check this my making sure that no SC atom is more than 20A (?) away from CA
    //needToRepack[i] = true;
    numeric::xyzVector< core::Real> ca_pos = fa_input_pose.residue(i).atom("CA").xyz();
    for ( Size j=1; j<=fa_input_pose.residue(i).natoms(); ++j ) {
      if ( (ca_pos - fa_input_pose.residue(i).atom(j).xyz()).length() > 20 ) {
	needToRepack[ i ] = true;
      }
    }

    //copy sidechains only for non-loop regions
    if ( !needToRepack[ i ] ) {
      pose.replace_residue(i, fa_input_pose.residue(i), true );
    }
  } //for 1..total_residue
} //copy_side_chains


bool ResolutionSwitcher::apply( core::pose::Pose &pose ) {
  bool const bCopySideChains( start_centroid_ || apply_to_centroid_ ); //we assume that start and apply pose are fullatom that side-chains are fine!
  utility::vector1< bool > needToRepack( pose.total_residue() , false );
  if ( apply_to_centroid_ ) {
    //puts full-atom sidechains on loop regions
    core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
    pose.conformation().detect_bonds();//apl fix this !
  }

  //find residues that have moved
  for ( Size i=1; i<=pose.total_residue(); ++i ) {
    if ( init_pose_.phi( i ) != pose.phi( i ) || init_pose_.psi( i ) != pose.psi( i ) ) {
      for ( Size j = i; j <= pose.total_residue() && j <= i + 3; j++ ) needToRepack[ j ] = true;
      //saftey buffer for repacking: 3 residues on each side of moved stuff
      if ( i > 1 && !needToRepack[ i-1 ] ) {
	for ( Size j = i-1; j >= 1 && j >= i-3; j++ ) needToRepack[ j ] = true;
      }
    }
  }

  if ( bCopySideChains && init_fa_ ) {
    copy_side_chains( pose, needToRepack, init_pose_ );
  }

   // repack loop + missing-density residues
  core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose );
  taskstd->restrict_to_repacking();
  taskstd->or_include_current(true);
  taskstd->restrict_to_residues(needToRepack);
  protocols::moves::PackRotamersMover pack1( scorefxn_fa , taskstd );
  pack1.apply( pose );

  // quick SC minimization
  core::optimization::AtomTreeMinimizer mzr;
  core::optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone", 1e-5, true, false );
  core::kinematics::MoveMap mm;
  mm.set_bb( false );
  mm.set_chi( true );
  mzr.run( pose, mm, *scorefxn_fa, options );

}



}

