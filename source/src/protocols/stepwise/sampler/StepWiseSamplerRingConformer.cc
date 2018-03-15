// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerRingConformer.cc
/// @brief Generate rotamer for one torsion angle.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerRingConformer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/rings/RingConformer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static basic::Tracer TR( "protocols.stepwise.sampler.StepWiseSamplerRingConformer" );

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerRingConformer::StepWiseSamplerRingConformer():
	ring_num_( 0 ),
	moving_rsd_( 0 ),
	ring_conformers_()
{}

StepWiseSamplerRingConformer::StepWiseSamplerRingConformer(
	Size const ring_num,
	Size const moving_rsd,
	core::chemical::rings::RingConformerSetCOP ring_conformer_set
):
	ring_num_( ring_num ),
	moving_rsd_( moving_rsd )
{
	ring_conformers_ = ring_conformer_set->get_all_nondegenerate_conformers();
}

StepWiseSamplerRingConformer::~StepWiseSamplerRingConformer()= default;

core::Size
StepWiseSamplerRingConformer::size() const { return ring_conformers_.size(); }

void StepWiseSamplerRingConformer::apply( core::pose::Pose & pose, Size const i ) {
	// How to ensure that i is appropriate?
	pose.set_ring_conformation( moving_rsd_, ring_num_, ring_conformers_[ i ] );
}
///////////////////////////////////////////////////////////////////////////
std::string StepWiseSamplerRingConformer::get_name() const {
	std::stringstream ss;
	ss << "StepWiseSamplerRingConformer residue:" << moving_rsd_ << " ring:" << ring_num_ << " conformers:" << ring_conformers_.size();
	return ss.str();
}

} //sampler
} //stepwise
} //protocols
