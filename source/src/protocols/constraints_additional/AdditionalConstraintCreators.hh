// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/protocols/constraints_additional/AdditionalConstraintCreators.hh
/// @brief  Base class for ConstraintCreators for the Constraint load-time factory registration scheme
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_constraints_additional_AdditionalConstraintCreators_hh
#define INCLUDED_protocols_constraints_additional_AdditionalConstraintCreators_hh

// Project Headers
#include <core/scoring/constraints/ConstraintCreator.hh>

// C/C++ headers
#include <string>

namespace protocols {
namespace constraints_additional {


/// @brief Mover creator for the BindingSiteConstraint constraint
class BindingSiteConstraintCreator : public core::scoring::constraints::ConstraintCreator
{
public:
	BindingSiteConstraintCreator();
	~BindingSiteConstraintCreator() override;

	core::scoring::constraints::ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Mover creator for the SequenceCoupling1BDconstraint
class SequenceCoupling1BDConstraintCreator : public core::scoring::constraints::ConstraintCreator
{
public:
	SequenceCoupling1BDConstraintCreator();
	~SequenceCoupling1BDConstraintCreator() override;

	core::scoring::constraints::ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Mover creator for the SequenceCouplingconstraint
class SequenceCouplingConstraintCreator : public core::scoring::constraints::ConstraintCreator
{
public:
	SequenceCouplingConstraintCreator();
	~SequenceCouplingConstraintCreator() override;

	core::scoring::constraints::ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief This class can be used to replace the standard AtomPairConstraintsCreator; see BrokerMain.cc for an example
class NamedAtomPairConstraintCreator : public core::scoring::constraints::ConstraintCreator
{
public:
	NamedAtomPairConstraintCreator();
	~NamedAtomPairConstraintCreator() override;

	core::scoring::constraints::ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};


} //namespace constraints_additional
} //namespace protocols

#endif
