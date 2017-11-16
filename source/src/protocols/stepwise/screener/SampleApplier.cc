// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/SampleApplier.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/sampler/StepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <protocols/moves/CompositionMover.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "protocols.stepwise.screener.SampleApplier" );

using namespace core;
using namespace protocols::stepwise::sampler;
using namespace protocols::stepwise::sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

//constructor -- this is not in use, but has to be filled in.
SampleApplier::SampleApplier():
	pose_( *( new pose::Pose ) )
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
SampleApplier::get_update( sampler::StepWiseSamplerOP sampler ){
	if ( !apply_residue_alternative_sampler_ &&
			sampler->type() == toolbox::RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ) {
		RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
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
