// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/RotamerSetOperations/RigidBodyMoveRotSetOps.cc
/// @brief  classes for rigid body movement during rotamer packing
/// @author Florian Richter, floric@u.washington.edu, sep 2009

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/RigidBodyMoveRotSetOps.hh>

//Project headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pose/Pose.hh>
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

static thread_local basic::Tracer tr( "protocols.toolbox.RotamerSetOperations.RigidBodyMoveRotSetOps" );

void
RigidBodyMoveBaseRSO::alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::pack::task::PackerTask const & ptask,
	core::graph::GraphCOP packer_neighbor_graph,
	core::pack::rotamer_set::RotamerSet & rotamer_set
){
	using namespace core::pack::rotamer_set;

	if (rotamer_set.get_n_residue_types() != 1)
	{
		tr.Debug << "alter_rotamer_set with multiple rotamer types in source set: " << rotamer_set.get_n_residue_types() << std::endl;
	}

	if (ptask.being_designed(rotamer_set.resid()))
	{
		tr.Error << "alter_rotamer_set called at designed position: " << rotamer_set.resid() << std::endl;
	}


	// Get sequence position for the rotamer set.
	core::Size sequence_position = rotamer_set.resid();

	// Generate list of alternate RB confs.
	utility::vector1< core::conformation::ResidueCOP > rigid_body_confs =
		get_rigid_body_confs(pose, ptask, sequence_position);

	tr.Debug << "At seqpos " << sequence_position << " generated " << rigid_body_confs.size() << " alternate rb conformations." << std::endl;

	if (rigid_body_confs.size() == 0 ) return;

	bump_selector_.reset();

	// Initialize residuetype
	core::chemical::ResidueTypeCOP concrete_residue;
	if( rotamer_set.num_rotamers() != 0 )
	{
		concrete_residue = (*rotamer_set.begin())->type().get_self_ptr();
	}
	else
	{
		concrete_residue = pose.residue( sequence_position ).type().get_self_ptr();
	}

	for (core::Size i = 1; i <= rigid_body_confs.size(); i++)
	{
		runtime_assert( rigid_body_confs[i]->name3() == pose.residue( sequence_position ).name3() );
	}

	// Check casting once
	bool cast_succesful( dynamic_cast< core::pack::rotamer_set::RotamerSet_ *  > (& rotamer_set) );

	//we need a couple of things analogous to RotamerSet_::build_rotamers_for_concrete
	//some code duplication for now, maybe this can be moved to a common function
	//1. Regenerate chi sampling information for the rotamer set
	utility::vector1< utility::vector1< core::Real > > extra_chi_steps( concrete_residue->nchi() );
	int nneighbs( pose.energies().tenA_neighbor_graph().get_node( sequence_position )->num_neighbors_counting_self() );
	bool buried( nneighbs >= int(ptask.residue_task( sequence_position ).extrachi_cutoff()) );
	if( cast_succesful ) {
		core::pack::rotamer_set::RotamerSet_ & rotset ( static_cast< core::pack::rotamer_set::RotamerSet_ & > (rotamer_set) );
		for ( Size ii = 1; ii <= concrete_residue->nchi(); ++ii ) {
			rotset.set_extra_samples( ptask, nneighbs, ii,	concrete_residue, extra_chi_steps[ ii ] );
		}
	}
	//RotamerSet_ duplicate lines over
	
	// Attempt to create rotamer library for the given type.
	core::pack::dunbrack::SingleResidueRotamerLibraryCAP rotlib =
		core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( *concrete_residue );

	if(!rotlib.expired())
	{
		tr.Debug << "At seqpos " << sequence_position << " retrieved rotamer library." << std::endl;
	}
	else
	{
		tr.Debug << "At seqpos " << sequence_position << " no rotamer library." << std::endl;
	}

	utility::vector1< core::conformation::ResidueOP > new_rots;

	for( core::Size i = 1; i <= rigid_body_confs.size(); ++i )
	{
		//let's make sure the residue used to create the rotamer
		//has the same connections as currently in the pose
		core::conformation::Residue existing_residue( *rigid_body_confs[i] );
		existing_residue.copy_residue_connections_from( pose.residue( sequence_position ) );
		utility::vector1< core::conformation::ResidueOP > suggested_rotamers_this_rbconf;

		// Generate full list of candidate rotamers at the rb conf if a rotamer library is available
		// otherwise just add the additional rb conf.
		core::pack::dunbrack::SingleResidueRotamerLibraryCOP rotlib_op = rotlib.lock();
		if( rotlib_op )
		{
			rotlib_op->fill_rotamer_vector( pose, sfxn, ptask, packer_neighbor_graph, concrete_residue, existing_residue, extra_chi_steps, buried, suggested_rotamers_this_rbconf);
		}
		else
		{
			//TODO fordas Loop through the initialized rotamer set instead of the pose orientation?
			core::conformation::ResidueOP currot( new core::conformation::Residue( pose.residue( sequence_position ) ) );
			currot->orient_onto_residue( existing_residue );
			suggested_rotamers_this_rbconf.push_back( currot );
		}

		// Prune rotamer list with bump check if needed, otherwise add all candidate rotamers
		for( core::Size j = 1; j<= suggested_rotamers_this_rbconf.size(); ++j )
		{
			core::conformation::ResidueOP new_rot = suggested_rotamers_this_rbconf[j];

			if( ptask.bump_check() && cast_succesful )
			{
				core::pack::rotamer_set::RotamerSet_ & rotset ( static_cast< core::pack::rotamer_set::RotamerSet_ & > (rotamer_set) );
				core::PackerEnergy bumpenergy = rotset.bump_check( new_rot, sfxn, pose, ptask, packer_neighbor_graph );
				BumpSelectorDecision decision =  bump_selector_.iterate_bump_selector( bumpenergy );
				switch ( decision ) {
					case KEEP_ROTAMER :
						new_rots.push_back( new_rot );
						break;
					case DELETE_PREVIOUS_ROTAMER :
						runtime_assert ( new_rots.size() > 0 );
						new_rots[ new_rots.size() ] = new_rot;
						break;
					case DELETE_ROTAMER : // do nothing
						break;
				}
			}
			else
			{
				new_rots.push_back( new_rot );
			}
		}// loop over suggested rotamers
	} //loop over rigid body confs

	//finally, add the new rotamers
	for( core::Size i = 1; i <= new_rots.size(); ++i)
	{
		rotamer_set.add_rotamer( *new_rots[i] );
	}
		
	tr.Debug <<
		"At seqpos " << sequence_position << " " <<
		new_rots.size() << " rotamers were added to rotamer set that now contains " <<
		rotamer_set.num_rotamers() << " rotamers." << std::endl;
}


core::Real
RigidBodyMoveBaseRSO::increase_packer_residue_radius(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskCOP task, //the_task
	core::Size residue_index)
{
	return determine_largest_nbr_atom_distance(
			pose.residue( residue_index ),
			get_rigid_body_confs(pose, *task, residue_index) );
}

core::Real
RigidBodyMoveBaseRSO::determine_largest_nbr_atom_distance(
	core::conformation::Residue const & target_res,
	utility::vector1< core::conformation::ResidueCOP > alternate_confs)
{
	if( alternate_confs.size() > 0 )
	{
		runtime_assert( target_res.name3() == alternate_confs[1]->name3() );
	}

	Size nbr_atom( target_res.nbr_atom() );
	core::PointPosition center_pos( target_res.xyz( nbr_atom ) );
	core::Real max_sq_dist(0.0);
	for( Size i = 1; i <= alternate_confs.size(); ++i){
		core::Vector dist_vect( center_pos - alternate_confs[i]->xyz( nbr_atom ) );
		core::Real sq_dist( dist_vect.length_squared() );
		if( sq_dist > max_sq_dist ) max_sq_dist = sq_dist;
	}

	return std::sqrt( max_sq_dist );
}

RigidBodyMoveRSO::RigidBodyMoveRSO( core::Size seqpos )
	: parent(),
		seqpos_(seqpos)
{
	rigid_body_confs_.clear();
}

RigidBodyMoveRSO::RigidBodyMoveRSO( RigidBodyMoveRSO const & other )
	: parent( other ),
		seqpos_(other.seqpos_),
		rigid_body_confs_(other.rigid_body_confs_)
{}

core::pack::rotamer_set::RotamerSetOperationOP
RigidBodyMoveRSO::clone() const{
	return core::pack::rotamer_set::RotamerSetOperationOP( new RigidBodyMoveRSO( *this ) );
}

utility::vector1< core::conformation::ResidueCOP >
RigidBodyMoveRSO::get_rigid_body_confs(
	core::pose::Pose const & /*pose*/,
	core::pack::task::PackerTask const & /*ptask*/,
	core::Size residue_index)
{
	runtime_assert(residue_index == seqpos_);

	return rigid_body_confs_;
}

void
RigidBodyMoveRSO::set_rigid_body_confs(
	utility::vector1< core::conformation::ResidueCOP > const & rigid_body_confs)
{
	rigid_body_confs_ = rigid_body_confs;
}

} //namespace protocols
} //namespace toolbox
} //namespace rotamer_set_operations


