// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/Mutator.cc
/// @brief  Implementation of the %Mutator class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/Mutator.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>

// numeric headers
#include <numeric/random/random_permutation.hh>



static basic::Tracer TR( "protocols.ligand_evolution.Mutator" );


namespace protocols {
namespace ligand_evolution {

Mutator::Mutator( FragmentLibrary const& library, utility::vector1< core::Real > const& weights, core::Real min_similarity, core::Real max_similarity )
:
	min_similarity_( min_similarity ),
	max_similarity_( max_similarity ),
	library_( library ),
	reaction_or_reagent_(weights )
{}

std::string const& Mutator::name() const {
	return name_;
}

utility::vector1< Individual > Mutator::apply( utility::vector1< Individual > const& parents, core::Size n_offspring ) const {

	if ( parents.empty() ) {
		TR.Warning << "Received no parents to mutate. Will return empty offspring list." << std::endl;
		return utility::vector1< Individual >();
	}

	if ( parents.size() > n_offspring ) {
		TR.Warning << "Got more parents than desired offspring. Not all of them will be mutated." << std::endl;
	}

	utility::vector1< Individual > local_parents( parents );
	numeric::random::random_permutation( local_parents.begin(), local_parents.end() );

	utility::vector1< Individual > offspring;

	core::Size next_parent = 1;
	while ( offspring.size() < n_offspring ) {
		offspring.push_back( mutate( local_parents.at( next_parent ) ) );
		next_parent++;
		if ( next_parent > local_parents.size() ) {
			next_parent = 1;
		}
	}

	return offspring;
}

Individual Mutator::mutate( Individual const& individual ) const {

	TR.Debug << "Start mutation for " << individual.identifier() << std::endl;
	// For most parts we are now working on the LigandIdentifier. This holds all information the FragmentLibrary needs
	LigandIdentifier new_identifier = individual.identifier();

	// Select which feature should be mutated. LigandIdentifier looks like [reaction_index, reagent1_index, reagent2_index, ... , reagentN_index]
	// I checked the WeightedSampler, it returns a number between [1,weights.size()] which translates to 1=reaction and 2=reagent
	core::Size feature = reaction_or_reagent_.random_sample();

	// Now that we know which feature should be mutated, we need to change the new identifier
	// If we wish to change one reagent, then it is easy. We ask the library to give us the list of similar ligands and pick one with a weighted selection
	// Reactions are just as easy. We set a new one and let the library search for two new most similar ligands
	if ( feature != 1 ) {
		core::Size reaction_positions = library_.reaction_positions( new_identifier[ 1 ] );
		feature = numeric::random::random_range( 1, static_cast< int >( reaction_positions ) ) + 1;
		TR.Debug << "Mutate reagent at position " << feature - 1 << std::endl;
		// a reagent should be changed                                  id of the current reagent, reaction of current ligand, position
		ReagentSimilarityList reagents( library_.get_similar_reagents( new_identifier.at( feature ), new_identifier.at( 1 ), feature - 1 ) );
		core::Size first_possible_reagent = 0;
		core::Size last_possible_reagent = reagents.size();
		// iterate over list and check which area fulfills the min and max similarity criteria
		for ( core::Size ii( 1 ); ii <= last_possible_reagent; ++ii ) {
			if ( first_possible_reagent == 0 && reagents.at( ii ).second <= max_similarity_ ) {
				first_possible_reagent = ii;
			}
			if ( reagents.at( ii ).second < min_similarity_ ) {
				last_possible_reagent = ii - 1;
			}
		}
		// check how many reagents can be considered
		TR.Debug << "First possible index: " << first_possible_reagent << ", last possible index: " << last_possible_reagent << std::endl;
		core::Size selected_reagent;
		if ( first_possible_reagent == 0 ) {
			// [ x x x x x x X |max| |min| ]
			// no other reagent is different enough to fulfill the max criteria, so select the most different e.g. last element in list
			selected_reagent = reagents.size();
		} else if ( first_possible_reagent >= last_possible_reagent ) {
			// [ x x x |max| |min| X x x x ]
			// the first reagent different enough to fulfill the max criteria is either the only one to fulfill the min criteria
			// or it already violates the min criteria, so it gets selected for being the most similar that is different enough
			selected_reagent = first_possible_reagent;
		} else {
			// [ x x |max| x X x x |min| x ]
			// there is are two or more possible reagents that fulfill the criteria so one of them is uniformly selected
			selected_reagent = numeric::random::random_range( first_possible_reagent, last_possible_reagent );
		}
		TR.Debug << "Selected reagent " << reagents.at( selected_reagent ).first << " at index " << selected_reagent << " with similarity " << reagents[ selected_reagent ].second <<
			". Range was [" << first_possible_reagent << "," << last_possible_reagent <<"]." << std::endl;
		new_identifier.at( feature ) = reagents.at( selected_reagent ).first;
	} else {
		// the reaction should be mutated, so we randomly select a new one that differs from the current and change both reagents to their most similar counterparts in the new reaction
		TR.Debug << "Changed reaction from " << new_identifier.at( 1 );
		new_identifier.at( 1 ) = library_.random_reaction( { new_identifier.at( 1 ) } );
		TR.Debug << " to " << new_identifier.at( 1 );
		core::Size reaction_positions = library_.reaction_positions( new_identifier[ 1 ] );
		for ( core::Size ii( 2 ); ii <= new_identifier.size(); ++ii ) {
			if ( ii - 1 <= reaction_positions ) {
				if ( new_identifier.at( ii ) == 0 ) {
					//this happens if we mutate from a smaller reaction into a larger one
					// simply select a random valid reagent and map it into reaction
					new_identifier.at( ii ) = numeric::random::random_range( 1, static_cast< int >( library_.reagents_size() ) );
				}
				ReagentSimilarityList reagents(
					library_.get_similar_reagents(new_identifier.at(ii), new_identifier.at(1), ii - 1));
				new_identifier.at(ii) = reagents.at(1).first;
				TR.Debug << ", similarity new reagent " << ii - 1 << ": " << reagents.at(1).second;
			} else {
				new_identifier[ ii ] = 0;
				TR.Debug << ", skip reagent " << ii - 1;
			}
		}
		TR.Debug << "." << std::endl;
	}

	TR.Debug << individual.identifier() << " mutated to " << new_identifier << std::endl;

	return Individual( new_identifier, { individual.id() }, "mutate" );
}

void Mutator::set_min_similarity( core::Real min_similarity ) {
	min_similarity_ = min_similarity;
}

void Mutator::set_max_similarity( core::Real max_similarity ) {
	max_similarity_ = max_similarity;
}

}
}
