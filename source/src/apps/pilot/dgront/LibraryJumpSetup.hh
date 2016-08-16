// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef LIBRARYJUMPSETUP_H_
#define LIBRARYJUMPSETUP_H_

// Unit Headers
#include <protocols/jumping/JumpSetup.fwd.hh>

// Package Headers
#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/JumpSample.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
//#include <core/scoring/constraints/ConstraintForest.fwd.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>


namespace protocols {
namespace jumping {

class LibraryJumpSetup : public BaseJumpSetup {
public:
	LibraryJumpSetup( Size , SecondaryStructureOP,core::pose::PoseOP,int ,int,int,int);
	LibraryJumpSetup( Size , SecondaryStructureOP ,core::pose::PoseOP ,int ,int *,int*,int*,int*);

  std::string type_name() const {
    return "LibraryJumpSetup";
  }
        virtual
        JumpSample
        create_jump_sample( ) const;

private:
        Size nres_;
        int nJumps_;
         SecondaryStructureOP ss_def_;
         core::pose::PoseOP native_pose_;
         int *iResid_;
         int *jResid_;
         int *orientation_;
         int *pleating_;
};
} }
#endif /* LIBRARYJUMPSETUP_H_ */
