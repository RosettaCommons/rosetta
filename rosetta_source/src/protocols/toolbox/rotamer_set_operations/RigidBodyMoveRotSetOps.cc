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
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

static basic::Tracer tr("protocols.toolbox.RotamerSetOperations.RigidBodyMoveRotSetOps");

RigidBodyMoveRSO::RigidBodyMoveRSO( core::Size seqpos )
	: parent(),
		seqpos_(seqpos)
{
	rigid_body_confs_.clear();
}

RigidBodyMoveRSO::RigidBodyMoveRSO( RigidBodyMoveRSO const & other )
	: parent( other ),
		seqpos_(other.seqpos_),
		rigid_body_confs_(other.rigid_body_confs_),
		bump_selector_(other.bump_selector_)
{}


RigidBodyMoveRSO::~RigidBodyMoveRSO(){}

core::pack::rotamer_set::RotamerSetOperationOP
RigidBodyMoveRSO::clone() const{
	return new RigidBodyMoveRSO( *this );
}


/// @brief fairly simple: all we need to do is iterate over
/// the rotamers in the rotamer_set, superimpose each of them
/// onto all the internally stored rigid body confs, and then
/// add the newly generated rotamers to the rotamer_set
/// @details we should also have some safety checks to make
/// sure all the rotamers in the set are of the same residue type,
/// as the internally stored ones, i.e. that no stupid user set this
/// position to designing in the task
void
RigidBodyMoveRSO::alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::pack::task::PackerTask const & ptask,
	core::graph::GraphCOP packer_neighbor_graph,
	core::pack::rotamer_set::RotamerSet & rotamer_set
){
	using namespace core::pack::rotamer_set;
	if (rigid_body_confs_.size() == 0 ) return;
	bump_selector_.reset();
	core::chemical::ResidueTypeCOP concrete_residue;
	if( rotamer_set.num_rotamers() != 0 ) concrete_residue = (*rotamer_set.begin())->type();
	else concrete_residue = pose.residue( seqpos_ ).type();
	runtime_assert( concrete_residue->is_ligand() );
	runtime_assert( rigid_body_confs_[1]->name3() == pose.residue( seqpos_ ).name3() );
	core::pack::dunbrack::SingleResidueRotamerLibraryCAP rotlib =
		core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( *concrete_residue );
	//only do this once
	bool cast_succesful( dynamic_cast< core::pack::rotamer_set::RotamerSet_ *  > (& rotamer_set) );

	//we need a couple of things analogous to RotamerSet_::build_rotamers_for_concrete
	//some code duplication for now, maybe this can be moved to a common function
	utility::vector1< utility::vector1< Real > > extra_chi_steps( concrete_residue->nchi() );
	int nneighbs( pose.energies().tenA_neighbor_graph().get_node( seqpos_ )->num_neighbors_counting_self() );
	bool buried( nneighbs >= int(ptask.residue_task( seqpos_).extrachi_cutoff()) );
	if( cast_succesful ) {
		core::pack::rotamer_set::RotamerSet_ & rotset ( static_cast< core::pack::rotamer_set::RotamerSet_ & > (rotamer_set) );
		for ( Size ii = 1; ii <= concrete_residue->nchi(); ++ii ) {
			rotset.set_extra_samples( ptask, nneighbs, ii,	concrete_residue, extra_chi_steps[ ii ] );
		}
	}
	//RotamerSet_ duplicate lines over

	utility::vector1< core::conformation::ResidueOP > new_rots;

	for( core::Size i = 1; i <= rigid_body_confs_.size(); ++i ){
		//let's make sure the residue used to create the rotamer
		//has the same connections as currently in the pose
		core::conformation::Residue existing_residue( *rigid_body_confs_[i] );
		existing_residue.copy_residue_connections_from( pose.residue( seqpos_ ) );
		utility::vector1< core::conformation::ResidueOP > suggested_rotamers_this_rbconf;

		if( rotlib ) rotlib->fill_rotamer_vector( pose, sfxn, ptask, packer_neighbor_graph, concrete_residue, existing_residue, extra_chi_steps, buried, suggested_rotamers_this_rbconf);
		else {
			core::conformation::ResidueOP currot( new core::conformation::Residue( pose.residue( seqpos_ ) ) );
			currot->orient_onto_residue( existing_residue );
			suggested_rotamers_this_rbconf.push_back( currot );
		}

		for( core::Size j = 1; j<= suggested_rotamers_this_rbconf.size(); ++j ){
			core::conformation::ResidueOP new_rot = suggested_rotamers_this_rbconf[j];

			if( ptask.bump_check() && cast_succesful ){
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
			else new_rots.push_back( new_rot );
		}// loop over suggested rotamers
	} //loop over rigid body confs

	//finally, add the new rotamers
	for( core::Size i = 1; i <= new_rots.size(); ++i) rotamer_set.add_rotamer( *new_rots[i] );
	tr << "At seqpos " << seqpos_ << " " << new_rots.size() << " rotamers were added to rotamer set that now contains " << rotamer_set.num_rotamers() << " rotamers." << std::endl;
}


RigidBodyMoveRSO::Real
RigidBodyMoveRSO::increase_packer_residue_radius(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskCOP //the_task
) const
{
	return determine_largest_nbr_atom_distance( pose.residue( seqpos_ ) );
}

void
RigidBodyMoveRSO::set_rigid_body_confs(
	utility::vector1< core::conformation::ResidueCOP > const & rigid_body_confs
)
{
	rigid_body_confs_ = rigid_body_confs;
}

RigidBodyMoveRSO::Real
RigidBodyMoveRSO::determine_largest_nbr_atom_distance(
	core::conformation::Residue const & target_res
) const
{
	if( rigid_body_confs_.size() > 0 ) runtime_assert( target_res.name3() == rigid_body_confs_[1]->name3() );
	Size nbr_atom( target_res.nbr_atom() );
	core::PointPosition center_pos( target_res.xyz( nbr_atom ) );
	Real max_sq_dist(0.0);
	for( Size i = 1; i <= rigid_body_confs_.size(); ++i){
		core::Vector dist_vect( center_pos - rigid_body_confs_[i]->xyz( nbr_atom ) );
		Real sq_dist( dist_vect.length_squared() );
		if( sq_dist > max_sq_dist ) max_sq_dist = sq_dist;
	}
	//std::cout << " largest nbr atom distance determined to be " << std::sqrt( max_sq_dist ) << std::endl;
	return std::sqrt( max_sq_dist );
}


} //namespace protocols
} //namespace toolbox
} //namespace rotamer_set_operations


