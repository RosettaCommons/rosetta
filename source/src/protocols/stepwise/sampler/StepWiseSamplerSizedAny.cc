// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSizedAny.cc
/// @brief Aggregate multiple samplers for modeler from any one of them.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerSizedAny.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.StepWiseSamplerSizedAny" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
StepWiseSamplerSizedAny::StepWiseSamplerSizedAny():
	StepWiseSamplerSized(),
	size_( 0 )
{}

StepWiseSamplerSizedAny::~StepWiseSamplerSizedAny(){}

///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSizedAny::init() {
	runtime_assert( !rotamer_list_.empty() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->init();
	}

	size_list_.clear();
	size_ = 0;

	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		Real const curr_size = rotamer_list_[i]->size();
		size_list_.push_back( (Size)curr_size );
		size_ = Size(size_ + curr_size);
		if ( curr_size == 0 ) TR << "Got a null rotamer sampler!" << std::endl;
	}

	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSizedAny::reset() {
	runtime_assert( is_init() );
	if ( random() ) {
		++( *this );
	} else {
		id_ = 1;
		curr_state_ = std::make_pair( 1, 1 );
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSizedAny::operator++() {
	runtime_assert( not_end() );
	if ( random() ) {
		id_ = numeric::random::rg().random_range( 1, size() );
		curr_state_ = id2state( id() );
	} else {
		++id_;
		++curr_state_.second;
		if ( curr_state_.second > size_list_[curr_state_.first] ) {
			curr_state_.second = 1;
			++curr_state_.first;
		}
	}
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSizedAny::apply( Pose & pose ) {
	runtime_assert( is_init() );
	rotamer_list_[curr_state_.first]->apply( pose, curr_state_.second );
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSizedAny::apply( core::pose::Pose & pose, Size const id ) {
	runtime_assert( is_init() );
	std::pair<Size, Size> const state = id2state( id );
	rotamer_list_[state.first]->apply( pose, state.second );
}
///////////////////////////////////////////////////////////////////////////
std::pair<Size, Size>
StepWiseSamplerSizedAny::id2state( core::Size const id ) const {
	runtime_assert( is_init() );
	runtime_assert( id <= size() );

	Size rotamer, id_rotamer( id );
	for ( rotamer = 1; rotamer < size_list_.size(); ++rotamer ) {
		if ( id_rotamer > size_list_[rotamer] ) {
			id_rotamer -= size_list_[rotamer];
		} else {
			break;
		}
	}

	return std::make_pair( rotamer, id_rotamer );
}
///////////////////////////////////////////////////////////////////////////
void StepWiseSamplerSizedAny::show( std::ostream & out, Size const indent ) const {
	StepWiseSamplerSized::show( out, indent );
	for ( Size n = 1; n <= indent+1; n++ ) out << ' ';
	out << "Following in multiple copies [" << rotamer_list_.size() << ']' << std::endl;
	rotamer_list_[1]->show( out, indent + 1 );
}
///////////////////////////////////////////////////////////////////////////
} //sampler
} //stepwise
} //protocols
