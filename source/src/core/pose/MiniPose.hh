// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/MiniPose.hh
/// @brief  MiniPose class
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_core_pose_MiniPose_HH
#define INCLUDED_core_pose_MiniPose_HH

// type headers
#include <core/pose/MiniPose.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

namespace core {
namespace pose {

/// @brief lightweight version of the pose with stuff I need.
// @details A very simple bag to hold xyz and fold_tree and sequence only. Kind of like a silent struct, but
/// not quite (silent_structs don't necessarily have xyz). Kind of like a Conformation but without the
/// overtree of an atom_tree.
/// @remark Should save memory compared to keeping the full pose (which includes atom_tree, energies etc.)
/// This is a bit like the SilentStruct -- although that class
/// has gotten a bit complicated -- easier to start from scratch.
class MiniPose : public utility::pointer::ReferenceCount  {

public:
	MiniPose( core::pose::Pose const & pose );

	MiniPose( utility::vector1< utility::vector1< PointPosition > > const & coords,
		core::kinematics::FoldTree const & fold_tree,
		std::string const & sequence );

	~MiniPose() override; // auto-removing definition from header{};

	core::kinematics::FoldTree const & fold_tree() const;

	utility::vector1< utility::vector1< PointPosition > > const & coords() const;

	utility::vector1< utility::vector1< std::string > > const & atom_names_list() const;

	utility::vector1< utility::vector1< std::string > > const & variant_types_list() const;

	Size size() const;

	std::string const & sequence() const;

	PointPosition const & xyz( core::id::AtomID atom_id ) const;

	std::string const & atom_name( core::id::AtomID atom_id ) const;

	utility::vector1< std::string > const & variant_types( Size const seq_num) const;

private: //data
	utility::vector1< utility::vector1< PointPosition > > coords_;
	core::kinematics::FoldTree fold_tree_;
	std::string sequence_;
	utility::vector1< utility::vector1< std::string > > atom_names_list_;
	utility::vector1< utility::vector1< std::string > > variant_types_list_;
};


} // namespace pose
} // namespace core

#endif // INCLUDED_core_pose_MiniPose_HH
