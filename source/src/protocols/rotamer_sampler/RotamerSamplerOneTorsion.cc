// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerOneTorsion.cc
/// @brief Generate rotamer for one torsion angle.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerOneTorsion.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerSamplerOneTorsion" );
static numeric::random::RandomGenerator RG( 2560199 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerSamplerOneTorsion::RotamerSamplerOneTorsion():
	RotamerSamplerOneValue(),
	torsion_id_()
{}

RotamerSamplerOneTorsion::RotamerSamplerOneTorsion(
		core::id::TorsionID const & tor_id,
		TorsionList const & allowed_torsions
):
	RotamerSamplerOneValue( allowed_torsions ),
	torsion_id_( tor_id )
{}

RotamerSamplerOneTorsion::~RotamerSamplerOneTorsion(){}

void RotamerSamplerOneTorsion::apply( core::pose::Pose & pose, Size const i ) {
	pose.set_torsion( torsion_id_, value( i ) );
}
///////////////////////////////////////////////////////////////////////////

} //rotamer_sampler
} //protocols
