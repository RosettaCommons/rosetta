// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/jump/JumpStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/jump/JumpStepWiseSampler.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility>

static basic::Tracer TR( "protocols.stepwise.sampler.jump.JumpStepWiseSampler" );

///////////////////////////////////////////////////////////////////
//
//  Bare bones JumpStepWiseSampler -- a more versatile version for general
//    rotation/translation searches, which
//    is also less sensitive to jump atoms, allows access to the 'Stub',
//    etc. is in rigid_body/RigidBodyStepWiseSampler.
//
//        -- rhiju, 2014
///////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampler {
namespace jump {

//Constructor
JumpStepWiseSampler::JumpStepWiseSampler():
	which_jump_( 0 )
{
	set_random( false );
}

//Constructor
JumpStepWiseSampler::JumpStepWiseSampler( Size const which_jump,
	utility::vector1< core::kinematics::Jump > const & jumps,
	bool const choose_random /* = false */ ):
	which_jump_( which_jump ),
	jumps_( jumps )
{
	set_random( choose_random );
}

//Destructor
JumpStepWiseSampler::~JumpStepWiseSampler() = default;

//////////////////////////////////////////////////////////////////////////
void
JumpStepWiseSampler::apply( core::pose::Pose & pose, Size const id )
{
	runtime_assert( id <= size() );
	pose.set_jump( which_jump_, jumps_[ id ] );
}


} //jump
} //sampler
} //stepwise
} //protocols
