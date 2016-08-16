// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/ScoringSchemeFactory.hh
/// @brief Factory for creating various types of constraints.
/// @author James Thompson <tex@u.washington.edu>

#ifndef INCLUDED_core_sequence_ScoringSchemeFactory_hh
#define INCLUDED_core_sequence_ScoringSchemeFactory_hh

#include <core/sequence/ScoringScheme.fwd.hh>
#include <map>

#include <iostream>


namespace core {
namespace sequence {

class ScoringSchemeFactory {
public:
	ScoringSchemeFactory(void);

	/// @brief adds a ScoringSchemeOP
	void add_type( ScoringSchemeOP new_scheme );
	ScoringSchemeOP get_scoring_scheme( std::string const & type ) const;

	typedef std::map< std::string, ScoringSchemeOP > ScoringSchemeTypes;
	ScoringSchemeTypes scheme_types_;
};

} // sequence
} // core

#endif
