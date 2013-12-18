// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/genetic_algorithm/Mutate1Randomizer.cc
/// @brief  Class to mutate only a single position at a time in any given entity
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/genetic_algorithm/Mutate1Randomizer.hh>

// Package headers
#include <protocols/genetic_algorithm/Entity.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

#include <numeric/random/random.fwd.hh>

#include <algorithm> // std::copy

#include <utility/exit.hh>


namespace protocols {
namespace genetic_algorithm {

Mutate1Randomizer::~Mutate1Randomizer() {}

void Mutate1Randomizer::mutate( protocols::genetic_algorithm::Entity & entity )
{
	if ( mutation_rate() == 1.0 ) {
		for ( Size ii = 1; ii <= entity.traits().size(); ++ii ) {
			core::Size const n_mutation_choices( choices()[ ii ].size() );
			Size new_element_ind = static_cast< core::Size >( numeric::random::uniform() * n_mutation_choices ) + 1;
			entity.set_entity_element( ii, choices()[ ii ][ new_element_ind ] );
		}
	} else {
		Size pos_to_mutate = static_cast< Size > ( entity.traits().size() * numeric::random::uniform() ) + 1;
		core::Size const n_mutation_choices( choices()[ pos_to_mutate ].size() );
		Size new_element_ind = static_cast< core::Size >( numeric::random::uniform() * n_mutation_choices ) + 1;
		entity.set_entity_element( pos_to_mutate, choices()[ pos_to_mutate ][ new_element_ind ] );
	}
}

} // namespace genetic_algorithm
} // namespace protocols

