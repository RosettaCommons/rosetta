// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/Ab_util.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef INCLUDED_protocols_antibody2_Ab_util_hh
#define INCLUDED_protocols_antibody2_Ab_util_hh


// Rosetta Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
// C++ Headers



///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody2 {

void simple_one_loop_fold_tree(
                                   core::pose::Pose & pose,
                                   loops::Loop const & loop
                                   );
    
    
    
    
void simple_fold_tree(
                          core::pose::Pose & pose_in,
                          core::Size jumppoint1,
                          core::Size cutpoint,
                          core::Size jumppoint2
                          );
    
    
void setup_simple_fold_tree(
                                core::Size jumppoint1,
                                core::Size cutpoint,
                                core::Size jumppoint2,
                                core::Size nres,
                                core::pose::Pose & pose_in );


} //namespace antibody2
} //namespace protocols


#endif //INCLUDED_protocols_loops_Ab_util_HH




