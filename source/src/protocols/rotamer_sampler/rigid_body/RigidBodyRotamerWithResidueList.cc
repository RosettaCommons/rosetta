// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueListRotamer.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.rigid_body.RigidBodyRotamerWithResidueList" );

using namespace protocols::rotamer_sampler::copy_dofs;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	//Constructor
	RigidBodyRotamerWithResidueList::RigidBodyRotamerWithResidueList( ResidueListRotamerOP copy_dofs_rotamer,
																																		RigidBodyRotamerOP rigid_body_rotamer ):
		copy_dofs_rotamer_( copy_dofs_rotamer ),
		rigid_body_rotamer_( rigid_body_rotamer )
	{
		// inner-most loop
		add_rotamer( copy_dofs_rotamer_ );
		// outer-most loop
		add_rotamer( rigid_body_rotamer_ );
	}

	//Destructor
	RigidBodyRotamerWithResidueList::~RigidBodyRotamerWithResidueList()
	{}

	/// @brief Apply the current rotamer to pose
	void
	RigidBodyRotamerWithResidueList::apply( core::pose::Pose & pose ){
		rigid_body_rotamer_->apply( pose, *copy_dofs_rotamer_->get_residue_at_origin() );
	}

	void
	RigidBodyRotamerWithResidueList::fast_forward_to_next_rigid_body(){
		copy_dofs_rotamer_->fast_forward();
	}

	void
	RigidBodyRotamerWithResidueList::fast_forward_to_next_translation(){
		copy_dofs_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_next_translation();
	}

	void
	RigidBodyRotamerWithResidueList::fast_forward_to_next_euler_gamma(){
		copy_dofs_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_next_euler_gamma();
	}

	void
	RigidBodyRotamerWithResidueList::fast_forward(){
		rigid_body_rotamer_->fast_forward_to_end();
		copy_dofs_rotamer_->fast_forward();
	}


	ValueList const &
	RigidBodyRotamerWithResidueList::get_rigid_body_values(){
		return rigid_body_rotamer_->get_values();
	}

	// from rigid body rotamer
	core::kinematics::Stub
	RigidBodyRotamerWithResidueList::get_stub(){
		return rigid_body_rotamer_->get_stub();
	}

	// from residue list rotamer
	core::conformation::ResidueOP
	RigidBodyRotamerWithResidueList::get_residue_at_origin(){
		return copy_dofs_rotamer_->get_residue_at_origin();
	}

	ResidueListRotamerOP
	RigidBodyRotamerWithResidueList::copy_dofs_rotamer(){ return copy_dofs_rotamer_; }

	RigidBodyRotamerOP
	RigidBodyRotamerWithResidueList::rigid_body_rotamer(){ return rigid_body_rotamer_; }


} //rigid_body
} //rotamer_sampler
} //protocols
