// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSized.cc
/// @brief Aggregate multiple samplers for sampling from any one of them.
/// @author Rhiju Das

// Unit headers
#include <protocols/rotamer_sampler/RotamerSized.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG( 2220192 );  // Magic number
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerSized" );

using namespace core;

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerSized::RotamerSized():
	RotamerBase(),
	id_( 0 )
{
	RotamerBase::set_random( false );
}

RotamerSized::~RotamerSized(){}

///////////////////////////////////////////////////////////////////////////
void RotamerSized::init() {
	id_ = 0;
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSized::reset() {
	runtime_assert( is_init() );
	if ( random() && not_end() ) {
		++( *this );
	} else {
		id_ = 1;
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerSized::operator++() {
	runtime_assert( not_end() );
	if ( random() ) {
		id_ = RG.random_range( 1, size() );
	} else {
		++id_;
	}
}
///////////////////////////////////////////////////////////////////////////
bool RotamerSized::not_end() const {
	runtime_assert( is_init() );
	return ( size() > 0 && id_ <= size() );
}

}
}
