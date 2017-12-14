// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/StubApplier.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/StubApplier.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.StubApplier" );

using namespace core;
using namespace protocols::stepwise::sampler;
using namespace protocols::stepwise::sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
StubApplier::StubApplier( kinematics::Stub & stub ):
	stub_( stub )
{}

//Destructor
StubApplier::~StubApplier() = default;

///////////////////////////////////////////////////////////////////
void
StubApplier::get_update( sampler::StepWiseSamplerOP sampler ){

	if ( sampler->type() == toolbox::RIGID_BODY_WITH_RESIDUE_LIST ) {
		RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
		stub_ = rigid_body_rotamer_with_copy_dofs.get_stub();
		return;
	} else if ( sampler->type() == toolbox::RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ) {
		RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
		stub_ = rigid_body_rotamer_with_residue_alternatives.get_stub();
		return;
	}

	runtime_assert( sampler->type() == toolbox::RIGID_BODY );
	RigidBodyStepWiseSampler & rigid_body_rotamer = *( static_cast< RigidBodyStepWiseSampler * >( sampler.get() ) );
	stub_ = rigid_body_rotamer.get_stub();

}

} //screener
} //stepwise
} //protocols
