// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerOneTorsion.cc
/// @brief Generate rotamer for one torsion angle.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerOneTorsion.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.StepWiseSamplerOneTorsion" );

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerOneTorsion::StepWiseSamplerOneTorsion():
	StepWiseSamplerOneValue(),
	torsion_id_()
{}

StepWiseSamplerOneTorsion::StepWiseSamplerOneTorsion(
	core::id::TorsionID const & tor_id,
	TorsionList const & allowed_torsions
):
	StepWiseSamplerOneValue( allowed_torsions ),
	torsion_id_( tor_id )
{}

StepWiseSamplerOneTorsion::~StepWiseSamplerOneTorsion(){}

void StepWiseSamplerOneTorsion::apply( core::pose::Pose & pose, Size const i ) {
	pose.set_torsion( torsion_id_, value( i ) );
}
///////////////////////////////////////////////////////////////////////////
std::string StepWiseSamplerOneTorsion::get_name() const {
	return "StepWiseSamplerOneTorsion residue:" + utility::to_string( torsion_id_.rsd() )+" torsion:"+
		to_string(torsion_id_.type())+':'+
		utility::to_string( torsion_id_.torsion() );
}

} //sampler
} //stepwise
} //protocols
