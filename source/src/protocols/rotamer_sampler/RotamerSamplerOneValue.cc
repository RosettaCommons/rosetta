// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerOneValue.hh
/// @brief Base class for RotamerSamplerOneValue
/// @author  Rhiju Das (rhiju@stanford.edu)

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerOneValue.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerSamplerOneValue" );
static numeric::random::RandomGenerator RG( 2560188 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerSamplerOneValue::RotamerSamplerOneValue():
	RotamerSamplerSized()
{}

RotamerSamplerOneValue::RotamerSamplerOneValue(
		ValueList const & allowed_values
):
	RotamerSamplerSized(),
	values_( allowed_values )
{}

RotamerSamplerOneValue::~RotamerSamplerOneValue(){}

} //rotamer_sampler
} //protocols
