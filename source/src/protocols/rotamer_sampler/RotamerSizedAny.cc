// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSizedAny.cc
/// @brief Aggregate multiple samplers for sampling from any one of them.
/// @detailed
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/RotamerSizedAny.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG( 256019 );  // Magic number
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerSizedAny" );

using namespace core;

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerSizedAny::RotamerSizedAny():
	RotamerSized(),
	size_( 0 ),
	id_( 0 )
{}

RotamerSizedAny::RotamerSizedAny( RotamerSizedAny const & other ):
	RotamerSized( other ),
	size_( other.size_ ),
	id_( other.id_ ),
	curr_state_( other.curr_state_ ),
	size_list_( other.size_list_ ),
	rotamer_list_( other.rotamer_list_ )
{}

RotamerSizedAny& RotamerSizedAny::operator=(
	RotamerSizedAny const &rhs
) {
	if ( this == &rhs ) return *this;
	RotamerSized::operator=( rhs );
	size_ = rhs.size_;
	id_ = rhs.id_;
	curr_state_ = rhs.curr_state_;
	size_list_ = rhs.size_list_;
	rotamer_list_ = rhs.rotamer_list_;
	return *this;
}

RotamerSizedAny::~RotamerSizedAny(){}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedAny::init() {
	runtime_assert( !rotamer_list_.empty() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->init();
	}

	size_list_.clear();
	size_ = 0;

	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		Real const curr_size = rotamer_list_[i]->size();
		size_list_.push_back( curr_size );
		size_ += curr_size;
	}

	runtime_assert( size_ != 0 );
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedAny::reset() {
	runtime_assert( is_init() );
	if ( is_random() ) {
		++( *this );
	} else {
		id_ = 1;
		curr_state_ = std::make_pair( 1, 1 );
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedAny::operator++() {
	runtime_assert( not_end() );
	if ( is_random() ) {
		id_ = RG.random_range( 1, size() );
		curr_state_ = id2state( id_ );
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
bool RotamerSizedAny::not_end() const {
	runtime_assert( is_init() );
	if ( is_random() ) return true;
	return id_ <= size();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedAny::apply( Pose & pose ) {
	runtime_assert( is_init() );
	rotamer_list_[curr_state_.first]->apply( pose, curr_state_.second );
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedAny::apply( core::pose::Pose & pose, Size const id ) {
	runtime_assert( is_init() );
	std::pair<Size, Size> const state = id2state( id );
	rotamer_list_[state.first]->apply( pose, state.second );
}
///////////////////////////////////////////////////////////////////////////
std::pair<Size, Size>
RotamerSizedAny::id2state( core::Size const id ) const {
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
}
}
