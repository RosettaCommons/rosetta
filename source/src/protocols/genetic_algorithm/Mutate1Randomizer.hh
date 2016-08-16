// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/genetic_algorithm/Mutate1Randomizer.hh
/// @brief  Class to mutate only a single position at a time in any given entity
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_genetic_algorithm_Mutate1Randomizer_hh
#define INCLUDED_protocols_genetic_algorithm_Mutate1Randomizer_hh

// Unit headers
#include <protocols/genetic_algorithm/Mutate1Randomizer.fwd.hh>

// Package headers
#include <protocols/genetic_algorithm/EntityRandomizer.hh>
#include <protocols/genetic_algorithm/Entity.fwd.hh>


namespace protocols {
namespace genetic_algorithm {

class Mutate1Randomizer : public PositionSpecificRandomizer {
public:
	typedef PositionSpecificRandomizer parent;

public:

	virtual ~Mutate1Randomizer();
	virtual void mutate( protocols::genetic_algorithm::Entity & entity );

};


} // namespace genetic_algorithm
} // namespace protocols

#endif
