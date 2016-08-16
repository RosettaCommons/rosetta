// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueListStepWiseSampler.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.rigid_body.RigidBodyStepWiseSamplerWithResidueList" );

using namespace protocols::stepwise::sampler::copy_dofs;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

//Constructor
RigidBodyStepWiseSamplerWithResidueList::RigidBodyStepWiseSamplerWithResidueList( ResidueListStepWiseSamplerOP copy_dofs_rotamer,
	RigidBodyStepWiseSamplerOP rigid_body_rotamer ):
	copy_dofs_rotamer_( copy_dofs_rotamer ),
	rigid_body_rotamer_( rigid_body_rotamer )
{
	// inner-most loop
	add_external_loop_rotamer( copy_dofs_rotamer_ );
	// outer-most loop
	add_external_loop_rotamer( rigid_body_rotamer_ );
}

//Destructor
RigidBodyStepWiseSamplerWithResidueList::~RigidBodyStepWiseSamplerWithResidueList()
{}

/// @brief Apply the current rotamer to pose
void
RigidBodyStepWiseSamplerWithResidueList::apply( core::pose::Pose & pose ){
	rigid_body_rotamer_->apply( pose, *copy_dofs_rotamer_->get_residue_at_origin() );
}

void
RigidBodyStepWiseSamplerWithResidueList::fast_forward_to_next_rigid_body(){
	copy_dofs_rotamer_->fast_forward();
}

void
RigidBodyStepWiseSamplerWithResidueList::fast_forward_to_next_translation(){
	copy_dofs_rotamer_->fast_forward();
	rigid_body_rotamer_->fast_forward_to_next_translation();
}

void
RigidBodyStepWiseSamplerWithResidueList::fast_forward_to_next_euler_gamma(){
	copy_dofs_rotamer_->fast_forward();
	rigid_body_rotamer_->fast_forward_to_next_euler_gamma();
}

void
RigidBodyStepWiseSamplerWithResidueList::fast_forward(){
	rigid_body_rotamer_->fast_forward_to_end();
	copy_dofs_rotamer_->fast_forward();
}


ValueList const &
RigidBodyStepWiseSamplerWithResidueList::get_rigid_body_values(){
	return rigid_body_rotamer_->get_values();
}

// from rigid body rotamer
core::kinematics::Stub
RigidBodyStepWiseSamplerWithResidueList::get_stub(){
	return rigid_body_rotamer_->get_stub();
}

// from residue list rotamer
core::conformation::ResidueOP
RigidBodyStepWiseSamplerWithResidueList::get_residue_at_origin(){
	return copy_dofs_rotamer_->get_residue_at_origin();
}

ResidueListStepWiseSamplerOP
RigidBodyStepWiseSamplerWithResidueList::copy_dofs_rotamer(){ return copy_dofs_rotamer_; }

RigidBodyStepWiseSamplerOP
RigidBodyStepWiseSamplerWithResidueList::rigid_body_rotamer(){ return rigid_body_rotamer_; }


} //rigid_body
} //sampler
} //stepwise
} //protocols
