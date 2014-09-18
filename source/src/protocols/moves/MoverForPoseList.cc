// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/MoverForPoseList.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/moves/MoverForPoseList.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

static thread_local basic::Tracer TR( "protocols.moves.MoverForPoseList" );

namespace protocols {
namespace moves {

	//Constructor
	MoverForPoseList::MoverForPoseList():
		Mover()
	{}

	//Destructor
	MoverForPoseList::~MoverForPoseList()
	{}

	// Yup, that's all this is for.
	void
	MoverForPoseList::apply( utility::vector1< core::pose::PoseOP > & pose_list ){
		for ( core::Size n = 1; n <= pose_list.size(); n++ ) {
			apply( *( pose_list[n] ) );
		}
	}

	// Could save time and memory by not copying into viewer_pose (see above)
	void
	MoverForPoseList::apply( utility::vector1< core::pose::PoseOP > & pose_list,
													 core::pose::Pose & viewer_pose ){
		utility::vector1< core::pose::PoseOP > output_pose_list;
		for ( core::Size n = 1; n <= pose_list.size(); n++ ){
			viewer_pose = ( *pose_list[n] ); //set viewer_pose;
			apply( viewer_pose );
			output_pose_list.push_back( viewer_pose.clone() );
		}
		pose_list = output_pose_list;
	}


} //moves
} //protocols
