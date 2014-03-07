// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/SampleApplier.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/rotamer_sampler/RotamerBase.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>
#include <protocols/moves/CompositionMover.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.SampleApplier" );

using namespace core;
using namespace protocols::rotamer_sampler;
using namespace protocols::rotamer_sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

	pose::Pose blank_pose;

	//constructor
	SampleApplier::SampleApplier():
		pose_( blank_pose )
	{
	}

	//constructor
	SampleApplier::SampleApplier( pose::Pose & pose,
																bool const apply_residue_alternative_sampler /* = true */ ):
		pose_( pose ),
		apply_residue_alternative_sampler_( apply_residue_alternative_sampler )
	{
	}

	//Destructor
	SampleApplier::~SampleApplier()
	{}

	void
	SampleApplier::get_update( rotamer_sampler::RotamerBaseOP sampler ){
		if ( !apply_residue_alternative_sampler_ &&
				 sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ){
			RigidBodyRotamerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyRotamerWithResidueAlternatives * >( sampler.get() ) );
			rigid_body_rotamer_with_residue_alternatives.apply_rigid_body_only( pose_ );
			return;
		}
		sampler->apply( pose_ );
	}

	void
	SampleApplier::apply_mover( moves::CompositionMoverOP mover, Size const i, Size const j ){
		mover->apply( pose_, i, j );
	}

} //screener
} //stepwise
} //protocols
