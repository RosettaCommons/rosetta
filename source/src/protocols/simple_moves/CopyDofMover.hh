// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/CopyDofMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_simple_moves_CopyDofMover_HH
#define INCLUDED_protocols_simple_moves_CopyDofMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/CopyDofMover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/copydofs/CopyDofs.fwd.hh>
#include <core/types.hh>
#include <map>

using namespace core;

namespace protocols {
namespace simple_moves {

class CopyDofMover: public protocols::moves::Mover {

public:

	//constructor
	CopyDofMover( pose::Pose const & template_pose, std::map< Size, Size > res_map );

	//destructor
	~CopyDofMover();

public:

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "CopyDofMover"; }

	pose::copydofs::CopyDofsInfo const & copy_dofs_info( pose::Pose const & pose ) const;

	void
	set_copy_dofs_info(
		pose::Pose const & pose,
		pose::copydofs::CopyDofsInfo const & copy_dofs_info );

private:

	std::string
	pose_string(
		pose::Pose const & pose,
		utility::vector1< Size > const & res_list ) const;

	std::string
	pose_string(
		pose::Pose const & pose
	) const;

	bool
	check_for_precomputed_copy_dofs_info( pose::Pose const & pose ) const;

private:

	pose::Pose const & template_pose_;
	pose::MiniPose const template_mini_pose_;
	std::map< Size, Size > res_map_;
	bool backbone_only_;
	bool side_chain_only_;
	bool ignore_virtual_;

	bool use_hash_;
	std::string pose_string_; // for precomputation. probably should hash to int.
	std::map< std::string, pose::copydofs::CopyDofsInfo > copy_dofs_info_;
};

} //simple_moves
} //protocols

#endif
