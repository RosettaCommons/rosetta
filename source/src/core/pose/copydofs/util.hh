// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/copydofs/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_copydofs_util_HH
#define INCLUDED_core_pose_copydofs_util_HH

#include <core/pose/copydofs/CopyDofs.fwd.hh>
#include <core/pose/copydofs/CopyDofsInfo.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>

namespace core {
namespace pose {
namespace copydofs {

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief A very useful function that copies degrees of freedom from one pose to another. res_map defines how to map residue numbers from the large pose to the smaller "scratch" pose.
/// @author rhiju, 2009.
void
copy_dofs(
	pose::Pose & pose,
	MiniPose const & scratch_pose,
	core::pose::ResMap const & res_map );

void
copy_dofs_match_atom_names( //Parin Sripakdeevong Dec 27, 2011.
	pose::Pose & pose,
	MiniPose const & chunk_pose,
	core::pose::ResMap const & res_map );

void
copy_dofs_match_atom_names(
	pose::Pose & pose,
	Pose const & scratch_pose);

void
copy_dofs(
	pose::Pose & pose,
	pose::Pose const & scratch_pose );

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
	pose::Pose & pose,
	pose::Pose const & scratch_pose,
	core::pose::ResMap const & res_map );

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs_match_atom_names(
	pose::Pose & pose,
	Pose const & scratch_pose,
	core::pose::ResMap const & res_map,
	bool const backbone_only   = false,
	bool const side_chain_only = false ,
	bool const ignore_virtual  = true );


////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
	pose::Pose & pose,
	Pose const & scratch_pose,
	std::map < id::AtomID , id::AtomID > const & atom_id_map);

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
	pose::Pose & pose,
	MiniPose const & scratch_pose,
	std::map < id::AtomID , id::AtomID > const & atom_id_map);

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
	pose::Pose & pose,
	MiniPose const & scratch_pose,
	std::map < id::AtomID , id::AtomID > const & atom_id_map );

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
	pose::Pose & pose,
	MiniPose const & scratch_pose,
	std::map < id::AtomID , id::AtomID > const & atom_id_map,
	std::map< id::AtomID, Size > const & atom_id_domain_map  );

///////////////////////////////////////////////////////////////////
void
setup_atom_id_map(
	std::map < core::id::AtomID , core::id::AtomID > & atom_id_map,
	core::pose::ResMap const & res_map,
	core::pose::Pose const & pose );

///////////////////////////////////////////////////////////////////
void
setup_atom_id_map_match_atom_names( //June 16, 2011 Parin Sripakdeevong
	std::map < core::id::AtomID , core::id::AtomID > & atom_id_map,
	ResMap const & res_map,
	core::pose::Pose const & pose,
	MiniPose const & chunk_pose );

///////////////////////////////////////////////////////////////////
void
setup_atom_id_map_match_atom_names(
	std::map < core::id::AtomID , core::id::AtomID > & atom_id_map,
	ResMap const & res_map,
	core::pose::Pose const & pose,
	core::pose::Pose const & reference_pose,
	bool const backbone_only   = false,
	bool const side_chain_only = false,
	bool const ignore_virtual  = true );

///////////////////////////////////////////////////////////////////
void
apply_dofs( pose::Pose & pose,
	CopyDofsInfo const & copy_dofs_info,
	core::Real const dof_tolerance = 1.0e-5 );

std::map< id::AtomID, Size >
blank_atom_id_domain_map( pose::Pose const & pose );

} //copydofs
} //pose
} //core

#endif
