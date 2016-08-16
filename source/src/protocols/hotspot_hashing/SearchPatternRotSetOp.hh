// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/hotspot_hashing/SearchPatternRotSetOp.hh
/// @brief  Creates rigid body variants from search pattern during repacking.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_protocols_hotspot_hashing_SearchPatternRotSetOp_hh
#define INCLUDED_protocols_hotspot_hashing_SearchPatternRotSetOp_hh

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/RigidBodyMoveRotSetOps.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

//Project headers
#include <core/kinematics/Stub.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/SearchPatternRotSetOp.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace hotspot_hashing {

class SearchPatternRotSetOp : public protocols::toolbox::rotamer_set_operations::RigidBodyMoveBaseRSO
{
public:
	typedef protocols::toolbox::rotamer_set_operations::RigidBodyMoveBaseRSO parent;

	SearchPatternRotSetOp( SearchPatternOP pattern );
	SearchPatternRotSetOp( SearchPatternRotSetOp const & other );

	virtual
	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	virtual
	utility::vector1< core::conformation::ResidueCOP > get_rigid_body_confs(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & ptask,
		core::Size residue_index);

	virtual
	core::Real increase_packer_residue_radius(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP,
		core::Size residue_index
	);

private:
	utility::vector1<core::kinematics::Stub> search_stubs_;

};

class AddSearchPatternRotSetOp : public core::pack::task::operation::TaskOperation
{
public:
	AddSearchPatternRotSetOp(core::Size target_residue, SearchPatternOP pattern);
	~AddSearchPatternRotSetOp();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	void parse_tag( TagCOP tag , DataMap & );

	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

private:
	core::Size target_residue_;
	SearchPatternOP pattern_;
};

}
}

#endif

