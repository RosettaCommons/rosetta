// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A wrapper for a very particular AmbiguousConstraint of MultiConstraints 
/// @author Andrew Watkins (amw579@stanford.edu, October 2016)

#ifndef INCLUDED_core_scoring_constraints_BasePairConstraint_fwd_hh
#define INCLUDED_core_scoring_constraints_BasePairConstraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {

class BasePairConstraint;

typedef utility::pointer::shared_ptr< BasePairConstraint > BasePairConstraintOP;
typedef utility::pointer::shared_ptr< BasePairConstraint const > BasePairConstraintCOP;

}
}
}

#endif

