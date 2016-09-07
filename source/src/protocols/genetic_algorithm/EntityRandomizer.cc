// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file EntityRandomizer.hh
/// @brief controls the alteration of the traits that define Entity
/// @author ashworth

// Unit headers
#include <protocols/genetic_algorithm/EntityRandomizer.hh>

#include <protocols/genetic_algorithm/Entity.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

#include <numeric/random/random.fwd.hh>

#include <algorithm> // std::copy

#include <utility/exit.hh>


namespace protocols {
namespace genetic_algorithm {

EntityRandomizer::EntityRandomizer() : utility::pointer::ReferenceCount(), entity_length_(0), mutation_rate_(1.), entity_template_(/* 0 */) {}
EntityRandomizer::~EntityRandomizer() = default;
EntityCOP EntityRandomizer::entity_template() const { return entity_template_; }
void EntityRandomizer::set_entity_template( EntityCOP entity) { entity_template_ = entity; }

void EntityRandomizer::set_entity_length( core::Size length ) { entity_length_ = length; }


////////////////////////////////////////////////////////////////////////////////////////////////////

DiscreteRandomizer::~DiscreteRandomizer()= default;
void DiscreteRandomizer::add_choice( EntityElementOP const & choice ) { choices_.push_back( choice ); }
void DiscreteRandomizer::set_choices( EntityElements const & choices ) { choices_ = choices; }
EntityElements const & DiscreteRandomizer::choices() const { return choices_; }

////////////////////////////////////////////////////////////////////////////////////////////////////
PositionSpecificRandomizer::~PositionSpecificRandomizer() = default;
utility::vector1< EntityElements > const &
PositionSpecificRandomizer::choices() const { return choices_; }

////////////////////////////////////////////////////////////////////////////////////////////////////
// method definitions

Entity::OP
EntityRandomizer::random_entity()
{
	runtime_assert( entity_length_ );
	EntityOP entity;
	if ( entity_template_ ) {
		entity = entity_template_->clone();
	} else {
		entity = EntityOP( new Entity );
	}
	entity->set_traits_size( entity_length_ );
	// temporarily increase mutation rate to 100%
	core::Real current_rate( mutation_rate_ );
	set_mutation_rate( 1. );
	mutate( *entity );
	// restore mutation rate
	set_mutation_rate( current_rate );

	return entity;
}

/// @brief randomly swap [1, N-1] traits between two entities
void
EntityRandomizer::crossover( Entity & entity1, Entity & entity2 )
{
	EntityElements traits1( entity1.traits() ), traits2( entity2.traits() );

	// track of the total number of traits that can potentially be swapped
	core::Size num_remaining(traits1.size());
	// pick number of traits to swap that is in the interval [1, N-1]
	core::Size num_to_swap(numeric::random::random_range(1, traits1.size()-1));

	for ( auto
			it1( traits1.begin() ), it2( traits2.begin() ),
			end1( traits1.end() ), end2( traits2.end() );
			(it1 != end1) && (it2 != end2); ++it1, ++it2 ) {
		core::Real probability_of_swap = static_cast<core::Real>(num_to_swap)/static_cast<core::Real>(num_remaining);
		if ( probability_of_swap >= numeric::random::uniform() ) {
			EntityElementOP temp = *it1;
			*it1 = *it2;
			*it2 = temp;
			--num_to_swap;
		}
		--num_remaining;
	}

	entity1.set_traits(traits1);
	entity2.set_traits(traits2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
DiscreteRandomizer::mutate( Entity & entity )
{
	EntityElements traits( entity.traits() );
	core::Size const size_choices( choices_.size() );
	for ( auto & trait : traits ) {
		if ( this->mutation_rate() < numeric::random::uniform() ) continue;
		trait = choices_[ static_cast< core::Size >( numeric::random::uniform() * size_choices ) + 1 ];
	}
	entity.set_traits(traits);
}

core::Size
DiscreteRandomizer::library_size() const
{
	core::Size size( choices_.size() );
	for ( core::Size i(1), e( this->entity_length() ); i < e; ++i ) size *= size;
	return size;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
PositionSpecificRandomizer::append_choices( EntityElements const & choices )
{
	choices_.push_back( choices );
	this->set_entity_length( choices_.size() );
}

void
PositionSpecificRandomizer::mutate( Entity & entity )
{
	EntityElements traits( entity.traits() );
	assert( traits.size() == choices_.size() );
	utility::vector1< EntityElements >::const_iterator choice_it( choices_.begin() );
	for ( auto
			it( traits.begin() ), end( traits.end() ); it != end;
			++it, ++choice_it ) {
		if ( this->mutation_rate() < numeric::random::uniform() ) continue;
		core::Size const size_choices( choice_it->size() );
		*it = (*choice_it)[ static_cast< core::Size >( numeric::random::uniform() * size_choices ) + 1 ];
	}
	entity.set_traits(traits);
}

core::Size
PositionSpecificRandomizer::library_size() const
{
	core::Size size( choices_.front().size() );
	for ( auto
			it( choices_.begin() ), end( choices_.end() ); it != end; ++it ) {
		if ( it == choices_.begin() ) continue;
		size *= it->size();
	}
	return size;
}

} // namespace genetic_algorithm
} // namespace protocols

