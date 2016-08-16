// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/MoverForPoseList.hh
/// @brief: In addition to applying to a pose, this mover can apply to a vector1 of poses.
/// @details: Each derived class should define its own apply() statement.
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_moves_MoverForPoseList_HH
#define INCLUDED_protocols_moves_MoverForPoseList_HH

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverForPoseList.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace moves {

class MoverForPoseList: public Mover {

public:

	//constructor
	MoverForPoseList();

	//destructor
	~MoverForPoseList();

public:

	using protocols::moves::Mover::apply;

	virtual void apply( Pose & ) = 0;

	// just apply to each member of the pose_list.
	virtual void apply( utility::vector1< core::pose::PoseOP > & pose_list );

	// goes through pose_list, copies into the viewer_pose (useful for graphics),
	// and then returns output pose_list.
	virtual void apply( utility::vector1< core::pose::PoseOP > & pose_list,
		core::pose::Pose & viewer_pose );

private:

};

} //moves
} //protocols

#endif
