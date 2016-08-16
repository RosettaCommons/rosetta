// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairFunction.cc
/// @brief  Count pair base class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

// STL Headers
#include <iostream>

// Utility Headers
#include <utility/exit.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

Real const CountPairFunction::cp_half( 0.2f );

bool
CountPairFunction::operator() (
	int const /*at1*/,
	int const /*at2*/,
	Real & /*weight*/,
	Size & /*path_dist*/
) const
{
	std::cerr << "Count pair ERROR: base class invocation of operator().";
	std::cerr <<  "Should have been overridden by derived class";
	utility_exit();
	return true;
}


} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core
