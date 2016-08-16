// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/SequenceProfileConstraintCreators.hh
/// @brief  Base class for ConstraintCreators for the Constraint load-time factory registration scheme
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_constraints_SequenceProfileConstraintCreator_HH
#define INCLUDED_core_scoring_constraints_SequenceProfileConstraintCreator_HH

// Project Headers
#include <core/scoring/constraints/ConstraintCreator.hh>

namespace core {
namespace scoring {
namespace constraints {

/// @brief Constraint creator for the SequenceProfileConstraint constraint
class SequenceProfileConstraintCreator : public core::scoring::constraints::ConstraintCreator
{
public:
	SequenceProfileConstraintCreator();
	virtual ~SequenceProfileConstraintCreator();

	virtual core::scoring::constraints::ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

} //namespace constraints
} //namespace scoring
} //namespace core

#endif // INCLUDED_core_scoring_constraints_SequenceProfileConstraintCreator_HH
