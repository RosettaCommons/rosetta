// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/selectors/TournamentSelector.cc
/// @brief  Implementation of the %TournamentSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/selectors/TournamentSelector.hh>

// utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/stream_util.hh>

static basic::Tracer TR( "protocols.ligand_evolution.TournamentSelector" );


namespace protocols {
namespace ligand_evolution {

TournamentSelector::TournamentSelector( core::Size tournament_size, core::Real probability )
:
	tournament_size_( tournament_size )
{
	set_probability( probability );
}

utility::vector1< Individual > TournamentSelector::apply( Population& population, core::Size size, bool remove ) const {

	if ( population.size() < size ) {
		TR.Error << "Can't create a subset of size " << size << " from a population of size " << population.size() << std::endl;
		utility_exit_with_message( "Population is to small for the desired subset size" );
	}

	// fill an array with possible indices to be selected
	utility::vector1< core::Size > indices;
	for ( core::Size ii( 1 ); ii <= population.size(); ++ii ) {
		indices.push_back( ii );
	}

	utility::vector1< core::Size > selected_indices;

	// now I need to run n_select tournaments always removing one index from indices and adding it to the selected indices
	while ( selected_indices.size() != size ) {
		// 1 - create a random subset from indices
		numeric::random::random_permutation( indices );
		core::Size actual_tournament_size = tournament_size_;
		if ( indices.size() < actual_tournament_size ) {
			actual_tournament_size = indices.size();
		}
		utility::vector1< core::Size > tournament_attendies( indices.begin(), indices.begin() + actual_tournament_size );

		// 2 - sort this subset. Since the indices expect a sorted population I can basically sort by indices.
		//     The individual at pos 1 has a lower score than the individual at pos 2, 3, 4 and so on
		std::sort( tournament_attendies.begin(), tournament_attendies.end() );

		// 3 - starting from the first index, select one with probability probability_ and extract it from indices and place it in selected_indices
		bool selected_winner = false;
		while ( !selected_winner ) {
			for ( core::Size index : tournament_attendies ) {
				if ( numeric::random::uniform() <= probability_ ) {
					selected_winner = true;
					selected_indices.push_back( index );
					indices.erase( std::remove( indices.begin(), indices.end(), index ), indices.end() );
					break;
				}
			}
		}
	}

	TR.Debug << "Final selection: " << selected_indices << std::endl;

	population.sort();
	utility::vector1< Individual > selected_individuals;
	if ( remove ) {
		selected_individuals = population.remove_individuals( selected_indices );
	} else {
		selected_individuals = population.individuals( selected_indices );
	}

	return selected_individuals;
}

std::string const& TournamentSelector::name() const {
	return name_;
}

void TournamentSelector::set_tournament_size( core::Size size ) {
	tournament_size_ = size;
}

void TournamentSelector::set_probability( core::Real probability ) {
	assert( probability > 0 && probability <= 1 );
	probability_ = probability;
}

}
}
