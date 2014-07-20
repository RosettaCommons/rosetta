// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamerSamplerWithResidueList.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerSamplerWithResidueList.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerSampler.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueListRotamerSampler.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.rigid_body.RigidBodyRotamerSamplerWithResidueList" );

using namespace protocols::rotamer_sampler::copy_dofs;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	//Constructor
	RigidBodyRotamerSamplerWithResidueList::RigidBodyRotamerSamplerWithResidueList( ResidueListRotamerSamplerOP copy_dofs_rotamer,
																																		RigidBodyRotamerSamplerOP rigid_body_rotamer ):
		copy_dofs_rotamer_( copy_dofs_rotamer ),
		rigid_body_rotamer_( rigid_body_rotamer )
	{
		// inner-most loop
		add_external_loop_rotamer( copy_dofs_rotamer_ );
		// outer-most loop
		add_external_loop_rotamer( rigid_body_rotamer_ );
	}

	//Destructor
	RigidBodyRotamerSamplerWithResidueList::~RigidBodyRotamerSamplerWithResidueList()
	{}

	/// @brief Apply the current rotamer to pose
	void
	RigidBodyRotamerSamplerWithResidueList::apply( core::pose::Pose & pose ){
		rigid_body_rotamer_->apply( pose, *copy_dofs_rotamer_->get_residue_at_origin() );
	}

	void
	RigidBodyRotamerSamplerWithResidueList::fast_forward_to_next_rigid_body(){
		copy_dofs_rotamer_->fast_forward();
	}

	void
	RigidBodyRotamerSamplerWithResidueList::fast_forward_to_next_translation(){
		copy_dofs_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_next_translation();
	}

	void
	RigidBodyRotamerSamplerWithResidueList::fast_forward_to_next_euler_gamma(){
		copy_dofs_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_next_euler_gamma();
	}

	void
	RigidBodyRotamerSamplerWithResidueList::fast_forward(){
		rigid_body_rotamer_->fast_forward_to_end();
		copy_dofs_rotamer_->fast_forward();
	}


	ValueList const &
	RigidBodyRotamerSamplerWithResidueList::get_rigid_body_values(){
		return rigid_body_rotamer_->get_values();
	}

	// from rigid body rotamer
	core::kinematics::Stub
	RigidBodyRotamerSamplerWithResidueList::get_stub(){
		return rigid_body_rotamer_->get_stub();
	}

	// from residue list rotamer
	core::conformation::ResidueOP
	RigidBodyRotamerSamplerWithResidueList::get_residue_at_origin(){
		return copy_dofs_rotamer_->get_residue_at_origin();
	}

	ResidueListRotamerSamplerOP
	RigidBodyRotamerSamplerWithResidueList::copy_dofs_rotamer(){ return copy_dofs_rotamer_; }

	RigidBodyRotamerSamplerOP
	RigidBodyRotamerSamplerWithResidueList::rigid_body_rotamer(){ return rigid_body_rotamer_; }


} //rigid_body
} //rotamer_sampler
} //protocols
