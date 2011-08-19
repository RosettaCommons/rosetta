// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  MultiResolutionProtocol
///   if your protocol wants to switch between full-atom and centroid representation derive it from this one
///   functionality to copy side-chains for fixed residues from an initial fa-pose to an  post-centroid fa-pose is provided
///
///
/// @author Oliver Lange
///


#ifndef INCLUDED_protocols_MultiResolutionProtocol_hh
#define INCLUDED_protocols_MultiResolutionProtocol_hh

// Unit Headers

// Package Headers
#include <core/pose/Pose.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <core/io/silent/silent.fwd.hh>
//#include <protocols/jumping/KinematicControl.hh>

// Project Headers
#include <protocols/moves/BoolMover.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {

class ResolutionSwitcher : public moves::BoolMover {
public:
  ResolutionSwitcher( core::pose::Pose const&, bool fullatom /*input pose*/, bool start_centroid, bool apply_to_centroid );

   //provide a full-atom pose from which side-chains are copied for unmoved torsions
  void set_side_chain_pose( core::pose::Pose& pose ) {
    fa_start_pose_( new core::pose::Pose( pose ) );
  }


  //@brief if input was full-atom but we started (start_pose) from centroid we will copy side-chains
  // repacks all residues that have been moved between start_pose and pose
  virtual bool apply( core::pose::Pose& );

  //@brief gives a starting pose ( with respect to setting in start_centroid )
  core::pose::Pose start_pose() const;

private:
  bool apply_to_centroid_;

  core::pose::Pose init_pose_; //if empty there is no fullatom info available
  bool init_fa_;

  bool start_centroid_;
};

}

#endif

