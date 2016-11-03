// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerComb.cc
/// @brief Aggregate of multiple rotamer samplers for modeler combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerComb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.StepWiseSamplerComb" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerComb::StepWiseSamplerComb():
	StepWiseSamplerBase(),
	is_empty_( false )
{}

StepWiseSamplerComb::~StepWiseSamplerComb(){}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerComb::init() {
	runtime_assert( !rotamer_list_.empty() );
	is_empty_ = false;
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->init();
		if ( !rotamer->not_end() ) {
			TR << "Got a null rotamer sampler!" << std::endl;
			is_empty_ = true;
		}
	}
	set_random( random() );
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerComb::reset() {
	runtime_assert( is_init() );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->reset();
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerComb::operator++() {
	runtime_assert( not_end() );

	if ( random() ) {
		for ( auto const & rotamer : rotamer_list_ ) {
			++( *rotamer );
		}
	} else {
		for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
			++( *rotamer_list_[i] );
			if ( rotamer_list_[i]->not_end() ) {
				break;
			} else {
				if ( i < rotamer_list_.size() ) rotamer_list_[i]->reset();
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////
bool StepWiseSamplerComb::not_end() const {
	runtime_assert( is_init() );
	if ( is_empty_ ) return false;
	if ( random() ) return true;
	for ( auto const & rotamer : rotamer_list_ ) {
		if ( !rotamer->not_end() ) return false;
	}
	return true;
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerComb::apply( Pose & pose ) {
	runtime_assert( is_init() );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->apply( pose );
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerComb::set_random( bool const setting ) {
	StepWiseSamplerBase::set_random( setting );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->set_random( setting );
		runtime_assert( rotamer->random() == setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerComb::show( std::ostream & out, Size const indent ) const {
	StepWiseSamplerBase::show( out, indent );
	// reverse direction so that 'inner loop' is last.
	for ( Size k = rotamer_list_.size(); k >= 1; k-- ) rotamer_list_[k]->show( out, indent + 1 );
}


///////////////////////////////////////////////////////////////////////////
} //sampler
} //stepwise
} //protocols
