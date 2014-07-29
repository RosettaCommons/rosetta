// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerAny.cc
/// @brief Aggregate multiple samplers for modeler from any one of them.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerAny.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static basic::Tracer TR( "protocols.sampler.StepWiseSamplerAny" );
static numeric::random::RandomGenerator RG( 2565849 );  // Magic number

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerAny::StepWiseSamplerAny():
	StepWiseSamplerBase(),
	curr_rotamer_( 1 ),
	is_weighted_( false ),
	is_empty_( true ),
	has_empty_( false )
{}

StepWiseSamplerAny::~StepWiseSamplerAny(){}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerAny::init() {
	runtime_assert( !rotamer_list_.empty() );
	is_empty_ = true;
	has_empty_ = false;
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->init();
		if ( !rotamer_list_[i]->not_end() ) {
			TR << "Got a null rotamer sampler!" << std::endl;
			has_empty_ = true;
		} else {
			is_empty_ = false;
		}
	}

	if ( !weights_.empty() ) {
		runtime_assert( weights_.size() == rotamer_list_.size() );
		is_weighted_ = true;
		Real sum( 0 );
		cdf_.clear();
		for ( Size i = 1; i <= weights_.size(); ++i ) {
			sum += weights_[i];
			cdf_.push_back( sum );
		}
		for ( Size i = 1; i <= cdf_.size(); ++i ) {
			cdf_[i] /= sum;
		}
		cdf_.back() = 1; // To prevent potential float-point error
	}

	set_random( random() );
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerAny::reset() {
	runtime_assert( is_init() );
	curr_rotamer_ = 1;
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->reset();
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerAny::operator++() {
	runtime_assert( not_end() );

	if ( random() ) {
		if ( is_weighted_ ) {
			Real const rand = RG.uniform();
			for ( Size i = 1; i <= cdf_.size(); ++i ) {
				if ( rand <= cdf_[i] ) {
					curr_rotamer_ = i;
					break;
				}
			}
		} else {
			curr_rotamer_ = RG.random_range( 1, rotamer_list_.size() );
		}
		if ( has_empty_ && !rotamer_list_[curr_rotamer_]->not_end() ) ++( *this );
		++( *rotamer_list_[curr_rotamer_] );
	} else {
		++( *rotamer_list_[curr_rotamer_] );
		while ( !rotamer_list_[curr_rotamer_]->not_end() ) {
			++curr_rotamer_;
		}
	}
}
///////////////////////////////////////////////////////////////////////////
bool StepWiseSamplerAny::not_end() const {
	runtime_assert( is_init() );
	if ( is_empty_ ) return false;
	if ( random() ) return true;
	if ( rotamer_list_[curr_rotamer_]->not_end() ) return true;
	return false;
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerAny::apply( Pose & pose ) {
	runtime_assert( is_init() );
	rotamer_list_[curr_rotamer_]->apply( pose );
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerAny::set_random( bool const setting ) {
	StepWiseSamplerBase::set_random( setting );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->set_random( setting );
		runtime_assert( rotamer_list_[i]->random() == setting );
	}
}
///////////////////////////////////////////////////////////////////////////
} //sampler
} //stepwise
} //protocols
