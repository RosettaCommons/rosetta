// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/util.hh
/// @brief  Utilities for ligand docking
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/util.hh>

#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/types.hh>
#include <core/pose/chains_util.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/stream_util.hh> // AUTO IWYU For operator<<

// Boost Headers

namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.util" );

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

/// @brief Move the center of specified multiple chains to the desired_centroid
void
move_ligand_to_desired_centroid(
	utility::vector1<std::string> const & chains,
	core::Vector const & desired_centroid,
	core::pose::Pose & pose
){
	TR.Debug << "Start move_ligand_to_desired_centroid using specified chains " << chains << std::endl;

	utility::vector1<core::Size> chain_ids;
	utility::vector1<core::Size> jump_ids;

	//obtain chain IDs and jump IDs from pose
	for ( std::string const & chain: chains ) {
		// Glycan ligands may have multiple chain ids (e.g. 4, 5)
		// for one single chain (e.g. X) due to the effects of branching
		// Therefore it is necessary to account for multiple chain ids
		// defining a single chain
		utility::vector1< core::Size > chain_ids_for_chain =
			core::pose::get_chain_ids_from_chain( chain, pose );
		for ( core::Size const & chain_id: chain_ids_for_chain ) {
			chain_ids.push_back( chain_id );
		}
		// Grab the Jump id corresponding to this chain (e.g. X)
		// If ligand were a branched glycan (e.g. chain X --> chain ids [4,5])
		// then it should be the first chain id (here 4) that has a Jump id
		// If the ligand is a small molecule with one chain (e.g. X)
		// and one chain id (e.g. 2), then chain id 2 should have a Jump id
		core::Size jump_id = core::pose::get_jump_id_from_chain_id(chain_ids_for_chain.front(), pose);
		jump_ids.push_back(jump_id);
		TR.Debug << "chain " << chain << " corresponds to chain ids " << chain_ids << " and Jump id " << jump_id << std::endl;
	}

	core::Vector const ligand_centroid = protocols::geometry::centroid_by_chains(pose, chain_ids);
	core::Vector const trans_vec = desired_centroid - ligand_centroid;
	core::Real const trans_len = trans_vec.length();

	if ( trans_len > 1e-3 ) { // otherwise we get NaNs
		protocols::rigid::RigidBodyTransMover mover(pose, jump_ids[1]);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.freeze();
		mover.apply(pose); //moving original chain

		//applying same mover to other chains
		for (
				core::Size counter = 2;
				counter <= jump_ids.size();
				++counter
				) {
			mover.rb_jump(jump_ids[counter]);
			mover.apply(pose);
		}
	}
	TR.Debug << "Finish move_ligand_to_desired_centroid using specified chains " << chains << std::endl;
}

/// @brief Move the center of the object(s) downstream of jump_id to the desired_centroid
void
move_ligand_to_desired_centroid(
	core::Size const jump_id,
	core::Vector const & desired_centroid,
	core::pose::Pose & pose
){
	TR.Debug << "Start move_ligand_to_desired_centroid by Jump id " << jump_id << std::endl;
	core::Vector const ligand_centroid = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
	core::Vector const trans_vec = desired_centroid - ligand_centroid;
	core::Real const trans_len = trans_vec.length();
	if ( trans_len > 1e-3 ) { // otherwise we get NaNs
		protocols::rigid::RigidBodyTransMover mover(pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}
	TR.Debug << "Finish move_ligand_to_desired_centroid by Jump id " << jump_id << std::endl;
}

/// @brief Move the neighbor atom of the specified multiple chains to the desired_position
void
move_ligand_neighbor_to_desired_position(
	utility::vector1<std::string> const & chains,
	core::Vector const & desired_position,
	core::pose::Pose & pose
){
	for ( std::string const & chain: chains ) {
		move_ligand_neighbor_to_desired_position(chain, desired_position,pose);
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
	if ( residues.size() == 0 ) {
		utility_exit_with_message("Can't find chain " + chain+ " in pose!");
	}
	if ( residues.size() != 1 ) {
		TR << "Multiple residues in chain " << chain << " only using first for positioning. " << std::endl;
	}
	core::conformation::ResidueCOP residue( residues[1] );
	core::Vector const ligand_position = residue->xyz( residue->nbr_atom() );

	core::Vector const trans_vec = desired_position - ligand_position;
	core::Real const trans_len = trans_vec.length();
	if ( trans_len > 1e-3 ) { // otherwise we get NaNs
		core::Size jump_id = core::pose::get_jump_id_from_chain(chain, pose);
		protocols::rigid::RigidBodyTransMover mover(pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}

}


} //namespace ligand_docking
} //namespace protocols
