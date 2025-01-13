// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/offspring_factory/IdentityFactory.cc
/// @brief  Implementation of the %IdentityFactory class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// unit headers
#include <protocols/ligand_evolution/offspring_factory/IdentityFactory.hh>

// utility headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.ligand_evolution.IdentityFactory" );


namespace protocols {
namespace ligand_evolution {

std::string const& IdentityFactory::name() const {
	return name_;
}

utility::vector1< Individual > IdentityFactory::apply( utility::vector1< Individual > const& parents, core::Size size ) const {

	if ( parents.size() < size ) {
		TR.Error << "Received " << parents.size() << " parents but should produce " << size << " offspring. Parent duplicates would be included." << std::endl;
        utility_exit_with_message("Copying individuals more than once leads to unexpected behavior. This should be avoided at all costs.");
	} else if ( parents.size() > size ) {
		TR.Warning << "Received " << parents.size() << " parents but should produce " << size << " offspring. Not all parents will be passed on." << std::endl;
	}

	utility::vector1< Individual > new_individuals;
	while ( new_individuals.size() < size ) {
		for ( core::Size ii( 1 ); ii <= parents.size() && new_individuals.size() < size; ++ii ) {
			new_individuals.emplace_back( parents.at( ii ) );
		}
	}

	return new_individuals;
}

}
}
