// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/NoOpStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/sampler/NoOpStepWiseSampler.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.NoOpStepWiseSampler" );

namespace protocols {
namespace stepwise {
namespace sampler {

//Constructor
NoOpStepWiseSampler::NoOpStepWiseSampler()
{
	init();
}

//Destructor
NoOpStepWiseSampler::~NoOpStepWiseSampler()
{}

} //sampler
} //stepwise
} //protocols
