// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneValue.hh
/// @brief Base class for StepWiseSamplerOneValue
/// @author  Rhiju Das (rhiju@stanford.edu)

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static thread_local basic::Tracer TR( "protocols.sampler.StepWiseSamplerOneValue" );

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerOneValue::StepWiseSamplerOneValue():
	StepWiseSamplerSized()
{}

StepWiseSamplerOneValue::StepWiseSamplerOneValue(
		ValueList const & allowed_values
):
	StepWiseSamplerSized(),
	values_( allowed_values )
{}

StepWiseSamplerOneValue::StepWiseSamplerOneValue(
		ValueList const & allowed_values,
		std::string const & tag
):
	StepWiseSamplerSized(),
	values_( allowed_values ),
	tag_( tag )
{}

StepWiseSamplerOneValue::~StepWiseSamplerOneValue(){}

std::string
StepWiseSamplerOneValue::get_name() const {
	std::string name = "StepWiseSamplerOneValue";
	if ( tag_.size() > 0 ) name += " tag:" + tag_;
	return name;
}

} //sampler
} //stepwise
} //protocols
