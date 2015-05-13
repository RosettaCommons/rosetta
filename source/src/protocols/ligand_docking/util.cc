// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/util.hh
/// @brief  Utilities for ligand docking
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/util.hh>

#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/pose/Pose.hh>

#include <core/types.hh>
#include <core/pose/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer TR( "protocols.ligand_docking.util" );

/// @brief Move the center of specified chain to the desired_centroid
void move_ligand_to_desired_centroid(
	std::string const & chain,
	core::Vector const & desired_centroid,
	core::pose::Pose & pose
) {
	core::Size jump_id = core::pose::get_jump_id_from_chain(chain, pose);
	// We probably want some error checking here to make sure that moving the jump
	// doesn't also move some other part of the protein due to FoldTree setup.
	move_ligand_to_desired_centroid(jump_id, desired_centroid, pose);
}

/// @brief Move the center of the object(s) downstream of jump_id to the desired_centroid
void
move_ligand_to_desired_centroid(
		core::Size const jump_id,
		core::Vector const & desired_centroid,
		core::pose::Pose & pose
){
	core::Vector const ligand_centroid = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
	core::Vector const trans_vec = desired_centroid - ligand_centroid;
	core::Real const trans_len = trans_vec.length();
	if (trans_len > 1e-3) { // otherwise we get NaNs
		protocols::rigid::RigidBodyTransMover mover(pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}
}

/// @brief Move the neighbor atom of the specified chain to the desired_position
void move_ligand_neighbor_to_desired_position(
	std::string const & chain,
	core::Vector const & desired_position,
	core::pose::Pose & pose
) {
	core::Size chain_id( core::pose::get_chain_id_from_chain(chain, pose) );
	core::conformation::ResidueCOPs residues( core::pose::get_chain_residues(pose, chain_id) );
	if( residues.size() == 0 ) {
		utility_exit_with_message("Can't find chain " + chain+ " in pose!");
	}
	if( residues.size() != 1 ) {
		TR << "Multiple residues in chain " << chain << " only using first for positioning. " << std::endl;
	}
	core::conformation::ResidueCOP residue( residues[1] );
	core::Vector const ligand_position = residue->xyz( residue->nbr_atom() );

	core::Vector const trans_vec = desired_position - ligand_position;
	core::Real const trans_len = trans_vec.length();
	if (trans_len > 1e-3) { // otherwise we get NaNs
		core::Size jump_id = core::pose::get_jump_id_from_chain(chain, pose);
		protocols::rigid::RigidBodyTransMover mover(pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}

}


} //namespace ligand_docking
} //namespace protocols
