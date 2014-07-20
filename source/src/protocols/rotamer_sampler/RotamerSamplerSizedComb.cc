// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerSizedComb.cc
/// @brief Aggregate of multiple rotamer samplers for sampling combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerSizedComb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG( 250199 );  // Magic number
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerSamplerSizedComb" );

using namespace core;

///////////////////////////////////////////////////////////////////////////
// RotamerSampler Combination is a lot like nesting loops, or like
//  writing numbers in decimal format.
//
// Due to a historical choice, however, rotamer samplers that are
//  later in the list are 'external' to rotamer samplers earlier in the list.
//
// This order seems counter-intuitive (cf. writing numbers -- incrementing changes
//  the *last* decimal place first), and may be worth fixing.
//
///////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerSamplerSizedComb::RotamerSamplerSizedComb():
	RotamerSamplerSized(),
	size_( 0 )
{}

RotamerSamplerSizedComb::RotamerSamplerSizedComb( RotamerSamplerSizedOP outer_loop_rotamer, RotamerSamplerSizedOP inner_loop_rotamer ):
	size_( 0 ) //will be fixed on init()
{
	rotamer_list_.push_back( inner_loop_rotamer );
	rotamer_list_.push_back( outer_loop_rotamer );
	init();
}

RotamerSamplerSizedComb::~RotamerSamplerSizedComb(){}

///////////////////////////////////////////////////////////////////////////
void RotamerSamplerSizedComb::init() {
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
}
///////////////////////////////////////////////////////////////////////////
void RotamerSamplerSizedComb::reset() {
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
void RotamerSamplerSizedComb::operator++() {
	runtime_assert( not_end() );
	if ( random() ) {
		id_ = RG.random_range( 1, size() );
	} else {
		++id_;
	}
	id_list_ = id2list( id_ );
	update_rotamer_ids();
}
///////////////////////////////////////////////////////////////////////////
void RotamerSamplerSizedComb::update_rotamer_ids() {
	for ( Size i = 1; i <= id_list_.size(); ++i ) {
		rotamer_list_[i]->set_id( id_list_[i] );
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerSamplerSizedComb::apply( Pose & pose ) {
	apply( pose, id_ );
}
///////////////////////////////////////////////////////////////////////////
void RotamerSamplerSizedComb::apply( core::pose::Pose & pose, Size const id ) {
	runtime_assert( is_init() );
	utility::vector1<Size> new_id_list = id2list( id );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->apply( pose, new_id_list[i] );
	}
}
///////////////////////////////////////////////////////////////////////////
utility::vector1<Size>
RotamerSamplerSizedComb::id2list( core::Size const id ) const {
	runtime_assert( is_init() );

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
RotamerSamplerSizedComb::list2id( utility::vector1<core::Size> const & id_list ) const {
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
RotamerSamplerSizedComb::fast_forward( Size const sampler_number ){
	runtime_assert( sampler_number <= rotamer_list_.size() );
	for ( Size n = 1; n <= sampler_number; n++ ){
		id_list_[ n ] = rotamer_list_[ n ]->size();
	}
	id_ = list2id( id_list_ );
}

/// @brief Set the random sampling state
void
RotamerSamplerSizedComb::set_random( bool const setting ){
	RotamerSamplerBase::set_random( setting );
	for ( Size n = 1; n <= rotamer_list_.size(); n++ )		rotamer_list_[ n ]->set_random( setting );
}

} //rotamer_sampler
} //protocols
