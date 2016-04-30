// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/InsertPoseIntoPoseMover.hh
/// @brief   Base class for graftmovers
/// @author  Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_grafting_InsertPoseIntoPoseMover_HH
#define INCLUDED_protocols_grafting_InsertPoseIntoPoseMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/grafting/simple_movers/InsertPoseIntoPoseMover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

/// @brief Insert a whole pose into another.  Loops, linkers, whaterver. No modeling here.  Wrapper to utility function insert_pose_into_pose.
/// @details Residues between start + end should be deleted before using this mover if needed.
///
class InsertPoseIntoPoseMover : public  protocols::moves::Mover {

public:

	InsertPoseIntoPoseMover(bool copy_pdbinfo = false);
	InsertPoseIntoPoseMover(core::pose::Pose const & src_pose, core::Size res_start, core::Size res_end, bool copy_pdbinfo = false);

	InsertPoseIntoPoseMover( InsertPoseIntoPoseMover const & src);

	virtual ~InsertPoseIntoPoseMover();

	virtual void
	apply(core::pose::Pose & pose);


public:
	void
	src_pose(core::pose::Pose const & src_pose);

	void
	start(core::Size res_start);

	core::Size
	start() const;

	void
	end(core::Size res_end);

	core::Size
	end() const;

public:
	virtual std::string
	get_name() const;

	protocols::moves::MoverOP
	clone() const;

	virtual void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	);

private:
	core::Size start_;
	core::Size end_;
	bool copy_pdbinfo_;
	core::pose::PoseOP src_pose_;
	TagCOP tag_; //This is so pdb_num can be parsed at apply time instead of construction time.
};


}
}
}


#endif  // INCLUDED_protocols_grafting_InsertPoseIntoPoseMover_HH
