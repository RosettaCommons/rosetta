// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneDOF.cc
/// @brief Generate rotamer for one DOF angle.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneDOF.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static basic::Tracer TR( "protocols.stepwise.sampler.StepWiseSamplerOneDOF" );

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerOneDOF::StepWiseSamplerOneDOF():
	StepWiseSamplerOneValue(),
	DOF_id_()
{}

StepWiseSamplerOneDOF::StepWiseSamplerOneDOF(
	core::id::DOF_ID const & tor_id,
	TorsionList const & allowed_DOFs
):
	StepWiseSamplerOneValue( allowed_DOFs ),
	DOF_id_( tor_id )
{}

StepWiseSamplerOneDOF::~StepWiseSamplerOneDOF()= default;

void StepWiseSamplerOneDOF::apply( core::pose::Pose & pose, Size const i ) {
	pose.set_dof( DOF_id_, value( i ) );
}
///////////////////////////////////////////////////////////////////////////
std::string StepWiseSamplerOneDOF::get_name() const {
	return "StepWiseSamplerOneDOF residue:" + utility::to_string( DOF_id_.rsd() )+" DOF:"+
		to_string(DOF_id_.type())+':'+
		utility::to_string( DOF_id_ );
}

} //sampler
} //stepwise
} //protocols
