// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef SPECIFICJUMPSETUP_H_
#define SPECIFICJUMPSETUP_H_

// Unit Headers
#include <protocols/jumping/JumpSetup.fwd.hh>

// Package Headers
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/JumpSample.hh>
#include <core/conformation/SecondaryStructure.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
//#include <core/scoring/constraints/ConstraintForest.hh>


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

class SpecificJumpSetup : public BaseJumpSetup {
public:
	SpecificJumpSetup( Size , core::conformation::SecondaryStructureOP,core::pose::PoseOP,int ,int,int,int);
	SpecificJumpSetup( Size , core::conformation::SecondaryStructureOP ,core::pose::PoseOP ,int ,int *,int*,int*,int*);

  std::string type_name() const {
    return "SpecificJumpSetup";
  }

        virtual
        JumpSample
        create_jump_sample( ) const;

        virtual
        JumpSample
        clean_jump_sample( JumpSample ) const {
	  std::cerr << "ERROR!!! unimplemented method in SpecificJumpSetup.cc " << std::endl;
	}

        /// @brief returns an ordered FragSet that is compatible with the JumpSample
        /// default: generate jumps from ss-library according to JumpSample
        // if the movemap allows sampling of the down-residue but not the up-residue:
        // include a jump with torsions only for the "down" residue
        // if the movemap allows neither sampling of up or down, don't include the jump
        virtual core::fragment::FragSetOP
        generate_jump_frags( JumpSample const&, core::kinematics::MoveMap const& mm) const;
private:
		Size nres_;
		int nJumps_;
		core::conformation::SecondaryStructureOP ss_def_;
		core::pose::PoseOP native_pose_;
		int *iResid_;
		int *jResid_;
		int *orientation_;
		int *pleating_;
};
} }
#endif /* SPECIFICJUMPSETUP_H_ */
