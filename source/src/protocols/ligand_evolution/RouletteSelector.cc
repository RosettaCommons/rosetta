// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/RouletteSelector.cc
/// @brief  Implementation of the %RouletteSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/RouletteSelector.hh>

// utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/stream_util.hh>

static basic::Tracer TR( "protocols.ligand_evolution.RouletteSelector" );


namespace protocols{
namespace ligand_evolution{

    utility::vector1< Individual > RouletteSelector::apply( Population& population, core::Size size, bool remove ) const {

        if( population.size() < size ) {
            TR.Error << "Can't create a subset of size " << size << " from a population of size " << population.size() << std::endl;
            utility_exit_with_message( "Population is to small for the desired subset size" );
        }

        // Sorting is save to be called frequent due to internal logic checks in population
        population.sort();

        // since scores can be positive and negative and I try to minimize instead of optimize I need to do two transformation steps
        // 1) score -= lowest_score to base all of them on 0
        // 2) score = highest_score - score to give the lowest score the highest weight
        // This transformation is not optimal since it changes the relative proportions, but I accept that
        core::Real lowest_score = 1000000;
        core::Real highest_score = -1000000;
        if( consider_positive_ ) {
            for( Individual const& individual : population.individuals ()) {
                core::Real score = individual.score ();
                if( score < lowest_score ) {
                    lowest_score = score;
                }
                if( score > highest_score ) {
                    highest_score = score;
                }
            }
            TR.Debug << "Found lowest " << lowest_score << " and highest score " << highest_score << std::endl;
        }

        numeric::random::WeightedReservoirSampler< core::Size > reservoir( size );

        TR.Debug << "Created weights [";

        for( core::Size ii( 1 ); ii <= population.size(); ++ii ) {
            core::Real weight;
            if( consider_positive_ ) {
                // using positive and negative weights together with the lowest value getting the highest weight
                weight = ( highest_score - lowest_score ) - ( population.individual( ii ).score() - lowest_score );
                if( weight == 0 ) {
                    // this ensures inclusion
                    weight += 0.00000000001;
                }
            } else {
                weight = -1.0 * population.individual( ii ).score();
            }
            reservoir.consider_sample( ii, weight );
            TR.Debug << weight << " (from " << population.individual( ii ).score() << " ), ";
        }

        TR.Debug << "]" << std::endl;

        utility::vector1< core::Size > selected_indices;
        reservoir.samples( &selected_indices );

        TR.Debug << name() << " selected " << selected_indices << std::endl;

        if( selected_indices.size() < size ) {
            TR.Warning << "Only " << selected_indices.size() << " individuals were selected but " << size << " are required. This can cause undesired behavior." << std::endl;
            if( !consider_positive_ ) {
                TR.Warning << "Positive scoring individuals are not considered. This is the likely source for the problem above." << std::endl;
            }
        }

        if( remove ) {
            return population.remove_individuals( selected_indices );
        } else {
            return population.individuals( selected_indices );
        }
    }

    std::string const& RouletteSelector::name() const {
        return name_;
    }

    void RouletteSelector::consider_positive( bool consider_positive ) {
        consider_positive_ = consider_positive;
    }

}
}
