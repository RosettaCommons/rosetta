// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSized.cc
/// @brief Aggregate multiple samplers for modeler from any one of them.
/// @author Rhiju Das

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static thread_local basic::Tracer TR( "protocols.sampler.StepWiseSamplerSized" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerSized::StepWiseSamplerSized():
	StepWiseSamplerBase(),
	id_( 0 )
{
	StepWiseSamplerBase::set_random( false );
}

StepWiseSamplerSized::~StepWiseSamplerSized(){}

///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSized::init() {
	id_ = 0;
	set_init( true );
	//	reset();
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSized::reset() {
	runtime_assert( is_init() );
	if ( random() && not_end() ) {
		++( *this );
	} else {
		id_ = 1;
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSized::operator++() {
	runtime_assert( not_end() );
	if ( random() ) {
		id_ = numeric::random::rg().random_range( 1, size() );
	} else {
		++id_;
	}
}
///////////////////////////////////////////////////////////////////////////
bool StepWiseSamplerSized::not_end() const {
	runtime_assert( is_init() );
	return ( size() > 0 && id_ <= size() );
}

} //sampler
} //stepwise
} //protocols
