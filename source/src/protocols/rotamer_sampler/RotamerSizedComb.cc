// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSizedComb.cc
/// @brief Aggregate of multiple rotamer samplers for sampling combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/RotamerSizedComb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG( 250199 );  // Magic number
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerSizedComb" );

using namespace core;

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerSizedComb::RotamerSizedComb():
	RotamerSized(),
	size_( 0 )
{}

RotamerSizedComb::~RotamerSizedComb(){}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedComb::init() {
	runtime_assert( !rotamer_list_.empty() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->init();
	}

	size_list_.clear();
	id_list_.clear();
	size_ = 1;

	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		Real const curr_size = rotamer_list_[i]->size();
		id_list_.push_back( 0 );
		size_list_.push_back( curr_size );
		size_ *= curr_size;
		if ( curr_size == 0 ) TR << "Got a null rotamer sampler!" << std::endl;
	}
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedComb::reset() {
	runtime_assert( is_init() );
	if ( random() ) {
		++( *this );
	} else {
		for ( Size i = 1; i <= id_list_.size(); ++i ) {
			id_list_[i] = 1;
		}
		id_ = 1;
	}
	update_rotamer_ids();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedComb::operator++() {
	runtime_assert( not_end() );

	if ( random() ) {
		id_ = RG.random_range( 1, size() );
		id_list_ = id2list( id_ );
	} else {
		++id_;
		++id_list_[1];
		for ( Size i = 1; i < id_list_.size(); ++i ) {
			if ( id_list_[i] > size_list_[i] ) {
				id_list_[i] = 1;
				++id_list_[i+1];
			}
		}
	}
	update_rotamer_ids();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedComb::update_rotamer_ids() {
	for ( Size i = 1; i <= id_list_.size(); ++i ) {
		rotamer_list_[i]->set_id( id_list_[i] );
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedComb::apply( Pose & pose ) {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->apply( pose, id_list_[i] );
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerSizedComb::apply( core::pose::Pose & pose, Size const id ) {
	runtime_assert( is_init() );

	utility::vector1<Size> new_id_list = id2list( id );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->apply( pose, new_id_list[i] );
	}
}
///////////////////////////////////////////////////////////////////////////
utility::vector1<Size>
RotamerSizedComb::id2list( core::Size const id ) const {
	runtime_assert( is_init() );
	runtime_assert( id <= size() );

	utility::vector1<Size> new_id_list ( id_list_.size(), 1 );

	new_id_list[1] = id;
	for ( Size i = 1; i < new_id_list.size(); ++i ) {
		if ( new_id_list[i] > size_list_[i] ) {
			Size quot = new_id_list[i] / size_list_[i];
			Size const rem = new_id_list[i] % size_list_[i];
			if ( rem == 0 ) {
				new_id_list[i] = size_list_[i];
				--quot;
			} else {
				new_id_list[i] = rem;
			}
			new_id_list[i+1] += quot;
		} else {
			break;
		}
	}
	return new_id_list;
}
///////////////////////////////////////////////////////////////////////////
core::Size
RotamerSizedComb::list2id( utility::vector1<core::Size> const & id_list ) const {
	runtime_assert( is_init() );

	Size id( 1 );
	Size block_size_( 1 );
	for ( Size i = 1; i <= id_list.size(); ++i ){
		id += block_size_ * ( id_list[ i ] - 1);
		block_size_ *= size_list_[ i ];
	}

	return id;
}

/// @brief Move sampler to end.
void
RotamerSizedComb::fast_forward( Size const sampler_number ){
	runtime_assert( sampler_number <= rotamer_list_.size() );
	for ( Size n = 1; n <= sampler_number; n++ ){
		id_list_[ n ] = rotamer_list_[ n ]->size();
	}
	id_ = list2id( id_list_ );
}

/// @brief Set the random sampling state
void
RotamerSizedComb::set_random( bool const setting ){
	RotamerBase::set_random( setting );
	for ( Size n = 1; n <= rotamer_list_.size(); n++ )		rotamer_list_[ n ]->set_random( setting );
}

}
}
