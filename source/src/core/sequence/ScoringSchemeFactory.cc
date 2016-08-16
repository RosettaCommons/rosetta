// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/sequence/ScoringSchemeFactory.cc
/// @brief Factory for creating various types of ScoringSchemes for use
/// in sequence alignment.
/// @author James Thompson

// Unit headers
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

// Package headers
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/DPScoringScheme.hh>
#include <core/sequence/L1ScoringScheme.hh>
#include <core/sequence/ProfSimScoringScheme.hh>
#include <core/sequence/CompassScoringScheme.hh>
#include <core/sequence/CompositeScoringScheme.hh>
#include <core/sequence/ChemicalShiftScoringScheme.hh>

#include <utility/exit.hh>

namespace core {
namespace sequence {

void
ScoringSchemeFactory::add_type( ScoringSchemeOP new_scheme ) {
	std::string const type_name( new_scheme->type() );
	scheme_types_[ type_name ] = new_scheme;
}

ScoringSchemeOP
ScoringSchemeFactory::get_scoring_scheme(
	std::string const & type
) const {
	ScoringSchemeTypes::const_iterator iter = scheme_types_.find( type );
	if ( iter != scheme_types_.end() ) {
		return iter->second->clone();
	} else {
		utility_exit_with_message(
			"ScoringSchemeFactory: unknown ScoringScheme type: " + type
		);
		return NULL;
	}
}

ScoringSchemeFactory::ScoringSchemeFactory(void) {
	// initialization of ScoringSchemes which this factory knows how to
	// instantiate
	add_type( ScoringSchemeOP( new SimpleScoringScheme() ) );
	add_type( ScoringSchemeOP( new MatrixScoringScheme() ) );
	add_type( ScoringSchemeOP( new L1ScoringScheme() ) );
	add_type( ScoringSchemeOP( new DPScoringScheme() ) );
	add_type( ScoringSchemeOP( new ProfSimScoringScheme() ) );
	add_type( ScoringSchemeOP( new CompassScoringScheme() ) );
	add_type( ScoringSchemeOP( new CompositeScoringScheme() ) );
	add_type( ScoringSchemeOP( new ChemicalShiftScoringScheme() ) );
}

} // namespace sequence
} // namespace core
