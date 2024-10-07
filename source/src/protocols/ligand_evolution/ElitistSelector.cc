// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/ElitistSelector.cc
/// @brief  Implementation of the %ElitistSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/ElitistSelector.hh>

// utility headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.ligand_evolution.ElitistSelector" );


namespace protocols{
namespace ligand_evolution{

    utility::vector1< Individual > ElitistSelector::apply( Population& population, core::Size size, bool remove ) const {
        if( population.size() < size ) {
            TR.Error << "Can't create a subset of size " << size << " from a population of size " << population.size() << std::endl;
            utility_exit_with_message( "Population is to small for the desired subset size" );
        }

        // Since the population is sorted we just need to create an array of ascending indices
        utility::vector1< core::Size > indices;
        for( core::Size ii( 1 ); ii <= size; ++ii ) {
            indices.push_back( ii );
        }

        population.sort();

        if( remove ) {
            return population.remove_individuals( indices );
        } else {
            return population.individuals( indices );
        }
    }

    std::string const& ElitistSelector::name() const {
        return name_;
    }

}
}
