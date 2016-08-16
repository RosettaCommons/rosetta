// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mutate1Randomizer.fwd.hh
/// @brief fowrad declaration of the class that controls the alteration of the traits that define an Entity
/// @author ashworth

#ifndef INCLUDED_protocols_genetic_algorithm_Mutate1Randomizer_FWD_HH
#define INCLUDED_protocols_genetic_algorithm_Mutate1Randomizer_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace genetic_algorithm {

class Mutate1Randomizer;

typedef utility::pointer::shared_ptr< Mutate1Randomizer > Mutate1RandomizerOP;
typedef utility::pointer::shared_ptr< Mutate1Randomizer const > Mutate1RandomizerCOP;

} // namespace genetic_algorithm
} // namespace protocols

#endif
