// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * SpecificJumpSetup.cpp
 *
 *  Created on: Nov 17, 2008
 *      Author: dgront
 */

#include "SpecificJumpSetup.hh"


// Unit Headers
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/util.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>

#include <core/fragment/FrameList.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/fragment/OrderedFragSet.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

//numeric headers
#include <numeric/random/random.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <fstream>


namespace protocols {
namespace jumping {

using namespace core;
using namespace fragment;

FragSetOP
SpecificJumpSetup::generate_jump_frags(
	JumpSample const& /*jumps*/,
	kinematics::MoveMap const& /*mm*/
) const {
        OrderedFragSetOP frags = new OrderedFragSet;

        for(int i=0;i<nJumps_;i++) {
          JumpingFrameOP our_jump = generate_jump_frame( iResid_[i],jResid_[i], true );
          our_jump->steal( *native_pose_ );
          frags->add( our_jump );
        }
        return frags;
}


SpecificJumpSetup::SpecificJumpSetup( Size nres, SecondaryStructureOP ss,core::pose::PoseOP native_one,int iResidue,int jResidue,int orientation,int pleating) {
        nres_ = nres;
        ss_def_ = ss;
        nJumps_ = 1;
        native_pose_ = native_one;
        iResid_ = new int[1];
        jResid_ = new int[1];
        orientation_ = new int[1];
        pleating_ = new int[1];
        iResid_[0] = iResidue;
        jResid_[0] = jResidue;
        orientation_[0] = orientation;
        pleating_[0] = pleating;
}

SpecificJumpSetup::SpecificJumpSetup( Size nres, SecondaryStructureOP ss,core::pose::PoseOP native_one,int nJumps,int *iResidue,int *jResidue,int *orientation,int *pleating) {
        nres_ = nres;
        ss_def_ = ss;
        nJumps_ = nJumps;
        native_pose_ = native_one;
        iResid_ = new int[nJumps];
        jResid_ = new int[nJumps];
        orientation_ = new int[nJumps];
        pleating_ = new int[nJumps];
        for(int i=0;i<nJumps;i++) {
          iResid_[i] = iResidue[i];
          jResid_[i] = jResidue[i];
          orientation_[i] = orientation[i];
          pleating_[i] = pleating[i];
        }
}

JumpSample
SpecificJumpSetup::create_jump_sample( ) const {
  PairingsList jumps;
  // ANTI = 1, INWARDS = 2
  for(int i=0;i<nJumps_;i++) {
    Pairing my_jump( iResid_[i],jResid_[i], orientation_[i], pleating_[i]); // 5,62,1,1
    jumps.push_back(my_jump);
  }

    // this one is not safe due to FArray,  generate random fold-tree  ( as used by SheetBuilder )
  return JumpSample ( nres_,  jumps, ss_def_->loop_fraction() );
}

} //jumping
} //protocols

