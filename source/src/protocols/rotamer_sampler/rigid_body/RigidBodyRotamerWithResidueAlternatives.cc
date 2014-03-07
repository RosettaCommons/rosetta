// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerComb.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.rigid_body.RigidBodyRotamerWithResidueAlternatives" );

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {


	//Constructor
	RigidBodyRotamerWithResidueAlternatives::RigidBodyRotamerWithResidueAlternatives( ResidueAlternativeRotamerCombOP residue_alternatives_rotamer,
																																										RigidBodyRotamerOP rigid_body_rotamer ):
		residue_alternatives_rotamer_( residue_alternatives_rotamer ),
		rigid_body_rotamer_( rigid_body_rotamer )
	{
		// inner-most loop
		add_rotamer( residue_alternatives_rotamer_ );
		// outer-most loop
		add_rotamer( rigid_body_rotamer_ );
	}

	//Destructor
	RigidBodyRotamerWithResidueAlternatives::~RigidBodyRotamerWithResidueAlternatives()
	{}

	void
	RigidBodyRotamerWithResidueAlternatives::apply_rigid_body_only( pose::Pose & pose ){
	 	rigid_body_rotamer_->apply( pose );
	}

	void
	RigidBodyRotamerWithResidueAlternatives::fast_forward_to_next_rigid_body(){
		residue_alternatives_rotamer_->fast_forward();
	}

	void
	RigidBodyRotamerWithResidueAlternatives::fast_forward_to_next_residue_pair( Size const i, Size const j){
		residue_alternatives_rotamer_->fast_forward_to_next_residue_pair( i, j );
	}

	void
	RigidBodyRotamerWithResidueAlternatives::fast_forward_to_next_residue( Size const i ){
		residue_alternatives_rotamer_->fast_forward_to_next_residue( i );
	}

	void
	RigidBodyRotamerWithResidueAlternatives::fast_forward_to_next_translation(){
		residue_alternatives_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_next_translation();
	}

	void
	RigidBodyRotamerWithResidueAlternatives::fast_forward_to_next_euler_gamma(){
		residue_alternatives_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_next_euler_gamma();
	}

	void
	RigidBodyRotamerWithResidueAlternatives::fast_forward(){
		residue_alternatives_rotamer_->fast_forward();
		rigid_body_rotamer_->fast_forward_to_end();
	}


	ValueList const &
	RigidBodyRotamerWithResidueAlternatives::get_rigid_body_values(){
		return rigid_body_rotamer_->get_values();
	}

	// from rigid body rotamer
	core::kinematics::Stub
	RigidBodyRotamerWithResidueAlternatives::get_stub(){
		return rigid_body_rotamer_->get_stub();
	}


	// used by ChainClosableScreener.
	conformation::ResidueCOP
	RigidBodyRotamerWithResidueAlternatives::get_residue( Size const seqpos ){
		transformed_residues[ seqpos ] = get_residue_at_origin( seqpos ).clone();
		rigid_body_rotamer_->apply( *transformed_residues[seqpos] ); // will apply rotation if in moving partition.
		return transformed_residues[seqpos];
	}

	// used by ChainClosableScreener.
	Vector
	RigidBodyRotamerWithResidueAlternatives::get_xyz( Size const seqpos, std::string const atom_name  ){
		//	transformed_residues[ seqpos ] = get_residue_at_origin( seqpos ).clone();
		Vector xyz = get_residue_at_origin( seqpos ).xyz( atom_name );
		rigid_body_rotamer_->apply( xyz, seqpos ); // will apply rotation if in moving partition.
		return xyz;
	}

	// from residue list rotamer
	core::conformation::Residue const &
	RigidBodyRotamerWithResidueAlternatives::get_residue_at_origin(){
		return get_residue_at_origin( rigid_body_rotamer_->moving_res() );
	}

	// from residue list rotamer
	core::conformation::Residue const &
	RigidBodyRotamerWithResidueAlternatives::get_residue_at_origin( Size const seqpos ){
		if ( residue_alternatives_rotamer_->has_resnum( seqpos ) ){
			return residue_alternatives_rotamer_->get_residue_at_origin( seqpos );
		} else {
			return rigid_body_rotamer_->pose_at_origin()->residue( seqpos );
		}
	}

	ResidueAlternativeRotamerCombOP
	RigidBodyRotamerWithResidueAlternatives::residue_alternatives_rotamer(){ return residue_alternatives_rotamer_; }

	RigidBodyRotamerOP
	RigidBodyRotamerWithResidueAlternatives::rigid_body_rotamer(){ return rigid_body_rotamer_; }


} //rigid_body
} //rotamer_sampler
} //protocols
