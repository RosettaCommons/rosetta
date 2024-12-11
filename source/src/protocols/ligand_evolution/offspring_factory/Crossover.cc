// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/offspring_factory/Crossover.cc
/// @brief  Implementation of the %Crossover class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)


// unit headers
#include <protocols/ligand_evolution/offspring_factory/Crossover.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/random/reservoir_sample.hh>

static basic::Tracer TR( "protocols.ligand_evolution.Crossover" );


namespace protocols {
namespace ligand_evolution {

Crossover::Crossover( FragmentLibrary const& library )
:
	library_( library )
{}

std::string const& Crossover::name() const {
	return name_;
}

utility::vector1< Individual > Crossover::apply( utility::vector1< Individual > const& parents, core::Size n_offspring ) const {

	utility::vector1< Individual > offspring;

	if ( parents.size() < 2 ) {

		if ( parents.empty() ) {
			TR.Warning << "Received no parents for crossover. Will return empty offspring list." << std::endl;
			return utility::vector1< Individual >();
		}

		TR.Warning << "At least two parents are needed for a crossover. This will return copies of the single parent only." << std::endl;
		while ( offspring.size() < n_offspring ) {
			offspring.push_back( Individual( parents.front().identifier(), { parents.front().id(), parents.front().id() }, "crossover" ) );
		}
		return offspring;
	}

	if ( parents.size() / 2 > n_offspring ) {
		TR.Warning << "Expected " << n_offspring << " offspring but got " << parents.size() << " parents. Not all of them will be used." << std::endl;
	}

    // todo turn it into a list of indices to manipulate
	utility::vector1< Individual > local_parents( parents );
	numeric::random::random_permutation( local_parents.begin(), local_parents.end() );
	core::Size current_index = 1;

	while ( offspring.size() < n_offspring ) {
		if ( current_index + 1 <= local_parents.size() ) {
			// If an odd number of parents is present, the last one won't participate in a crossover for this iteration, since two parents are always required
			offspring.push_back( cross( local_parents.at( current_index ), local_parents.at( current_index + 1 ) ) );
			current_index += 2;
		} else {
			// shuffle parents again to allow for different crossovers in the next iteration
			current_index = 1;
            numeric::random::random_permutation( local_parents.begin(), local_parents.end() );
		}
	}

	return offspring;
}

Individual Crossover::cross( Individual const& reaction_parent, Individual const& other_parent ) const {

	// TODO Implement the possibility to get reagent1 as a new reagent2 similarity based like A-B + C-D = A-C'

	// start with the identifier from the reaction parent
	LigandIdentifier new_identifier = reaction_parent.identifier();

	// next create a list of features that should be extracted from the other parent
	// this can be as many features as it wants but min one needs to come from the reaction giving parent and one from itself, otherwise this would be point mutation on reaction
	// randomly determine how many features should be selected
	core::Size reaction_positions = std::min( library_.reaction_positions( reaction_parent.identifier().at( 1 ) ), library_.reaction_positions( other_parent.identifier().at( 1 ) ) );
	core::Size n_features = numeric::random::random_range( 1, static_cast< int >( reaction_positions - 1 ) );

	// features accepts values added to it until it has n_features. Afterwards it replaces new features probabilistic so that all values have a chance of 1/N to end up in features
	numeric::random::ReservoirSampler< core::Size > features( n_features );
	for ( core::Size feature( 2 ); feature <= reaction_positions + 1; ++feature ) {
		features.add_value( feature );
	}

	// these features will now be used to replace entries from the reaction parent with one from the other parent
	for ( core::Size feature : features.values() ) {
		LigandIdentifier const& other_identifier = other_parent.identifier();
		if ( other_identifier.at( 1 ) == new_identifier.at( 1 ) ) {
			// both parents come from the same reaction
			new_identifier.at( feature ) = other_identifier.at( feature );
		} else {
			// both parents come from different reactions and reagents have to be mapped. Select the most similar reagent.
			new_identifier.at( feature ) = library_.get_similar_reagents( other_identifier.at( feature ), new_identifier.at( 1 ), feature - 1 ).front().first;
		}
	}

	TR.Debug << reaction_parent.identifier() << " and " << other_parent.identifier() << " made " << new_identifier << std::endl;

	return Individual( new_identifier, { reaction_parent.id(), other_parent.id() }, "crossover" );
}

}
}
