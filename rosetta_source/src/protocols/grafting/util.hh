// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/grafting/util.hh
/// @brief Header for grafting utility functions.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_grafting_util_hh
#define	INCLUDED_protocols_grafting_util_hh

#include <core/pose/Pose.hh>


namespace protocols {
namespace grafting {
using core::pose::Pose;
using core::Size;
    
    
///@brief Deletes a region of the pose. Starting from and including 'start' and 'end' residue.
void 
delete_region(Pose & pose, Size const start, Size const end);
        
///@brief Returns a region of the pose including start and end as a new pose. Has a simple foldtree.
Pose 
return_region(Pose & pose, Size const start, Size const end);
    
///@brief replaces residues from from_pose to to_pose into pose where insertion region is defined. Returns product as a new value.
Pose
replace_region(Pose const & from_pose, Pose const & to_pose, Size const from_pose_start_residue, Size const to_pose_start_residue, Size const insertion_length);

///@author Steven Lewis smlewi@gmail.com
///@brief inserts one pose into another pose, returning the product as a new value. 
///@details Nter->Cter. Coordinates and dihedrals of insert are unchanged.
///@details Begins insertion AFTER insert point.

Pose
insert_pose_into_pose(Pose const & scaffold_pose, Pose const & insert_pose, Size const insert_point, Size const insert_point_end);
    
}//namespace grafting
}//namespace protocols


#endif	//INCLUDED_protocols_grafting_util_hh

