// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/count_pair/CountPairFunction.hh
/// @brief  Count pair base class interface
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_types_hh
#define INCLUDED_core_scoring_etable_count_pair_types_hh

#include <core/types.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

enum  CPCrossoverBehavior {
	CP_CROSSOVER_3 = 1,
	CP_CROSSOVER_3FULL,
	CP_CROSSOVER_4,
	CP_CROSSOVER_34,
	nCPCrossoverBehaviors = CP_CROSSOVER_34 //keep this guy last
};

enum CPResidueConnectionType {
	CP_NO_BONDS = 1,
	CP_ONE_BOND,
	CP_MULTIPLE_BONDS_OR_PSEUDOBONDS,
	nCPResidueConnectionTypes = CP_MULTIPLE_BONDS_OR_PSEUDOBONDS // keep this guy last
};

// if an atom is separated by at least 5 bonds from any other atom, it is
// effectively at an infinite separation
Size const INFINITE_SEPARATION( 5 );

}
}
}
}

#endif
