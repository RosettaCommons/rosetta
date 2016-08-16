// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/StubDistanceScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/StubDistanceScreener.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueList.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.screener.StubDistanceScreener" );

using namespace protocols::stepwise::sampler;
using namespace protocols::stepwise::sampler::rigid_body;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
StubDistanceScreener::StubDistanceScreener(  core::kinematics::Stub & moving_res_base_stub,
	core::kinematics::Stub const & reference_stub,
	core::Real const max_distance_squared ):
	moving_res_base_stub_( moving_res_base_stub ),
	reference_stub_( reference_stub ),
	max_distance_squared_( max_distance_squared )
{
}

//Destructor
StubDistanceScreener::~StubDistanceScreener()
{}


//////////////////////////////////////////
bool
StubDistanceScreener::check_screen() {
	// TR << ( moving_res_base_stub_.v - reference_stub_.v ).length() << " " <<  std::sqrt( max_distance_squared_ ) << std::endl;;
	return ( ( moving_res_base_stub_.v - reference_stub_.v ).length_squared() <= max_distance_squared_ );
}

//////////////////////////////////////////
void
StubDistanceScreener::fast_forward( sampler::StepWiseSamplerBaseOP sampler ){
	if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_LIST ) {
		RigidBodyStepWiseSamplerWithResidueList & rigid_body_rotamer_with_copy_dofs = *( static_cast< RigidBodyStepWiseSamplerWithResidueList * >( sampler.get() ) );
		rigid_body_rotamer_with_copy_dofs.fast_forward_to_next_translation();
	} else if ( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES ) {
		RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );
		rigid_body_rotamer_with_residue_alternatives.fast_forward_to_next_translation();
	} else if ( sampler->type() == RIGID_BODY ) {
		RigidBodyStepWiseSampler & rigid_body_rotamer = *( static_cast< RigidBodyStepWiseSampler * >( sampler.get() ) );
		rigid_body_rotamer.fast_forward_to_next_translation();
	}
}


} //screener
} //stepwise
} //protocols
