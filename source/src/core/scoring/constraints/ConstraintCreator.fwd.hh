// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ConstraintCreator.hh
/// @brief  Base class for ConstraintCreators for the Constraint load-time factory registration scheme
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_constraints_ConstraintCreator_fwd_hh
#define INCLUDED_core_scoring_constraints_ConstraintCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {

/// @brief Abstract base class for a Constraint factory; the Creator class is responsible for
/// creating a particular Constraint class.
class ConstraintCreator;

typedef utility::pointer::shared_ptr< ConstraintCreator > ConstraintCreatorOP;
typedef utility::pointer::shared_ptr< ConstraintCreator const > ConstraintCreatorCOP;

} //namespace constraints
} //namespace scoring
} //namespace core

#endif
