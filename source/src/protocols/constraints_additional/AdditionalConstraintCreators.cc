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


// Unit Headers
#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>

/// Project Headers
#include <protocols/constraints_additional/SequenceCouplingConstraint.hh>
#include <protocols/constraints_additional/SequenceCoupling1BDConstraint.hh>
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>
#include <protocols/constraints_additional/BindingSiteConstraint.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace constraints_additional {


BindingSiteConstraintCreator::BindingSiteConstraintCreator() {}
BindingSiteConstraintCreator::~BindingSiteConstraintCreator() {}

core::scoring::constraints::ConstraintOP
BindingSiteConstraintCreator::create_constraint() const {
	return core::scoring::constraints::ConstraintOP( new BindingSiteConstraint );
}

std::string BindingSiteConstraintCreator::keyname() const
{
	return "BindingSite";
}

SequenceCoupling1BDConstraintCreator::SequenceCoupling1BDConstraintCreator() {}
SequenceCoupling1BDConstraintCreator::~SequenceCoupling1BDConstraintCreator() {}

core::scoring::constraints::ConstraintOP
SequenceCoupling1BDConstraintCreator::create_constraint() const
{
	return core::scoring::constraints::ConstraintOP( new SequenceCoupling1BDConstraint );
}

std::string
SequenceCoupling1BDConstraintCreator::keyname() const
{
	return "SequenceCoupling1BD";
}
SequenceCouplingConstraintCreator::SequenceCouplingConstraintCreator() {}
SequenceCouplingConstraintCreator::~SequenceCouplingConstraintCreator() {}

core::scoring::constraints::ConstraintOP
SequenceCouplingConstraintCreator::create_constraint() const
{
	return core::scoring::constraints::ConstraintOP( new SequenceCouplingConstraint );
}

std::string
SequenceCouplingConstraintCreator::keyname() const
{
	return "SequenceCoupling";
}

NamedAtomPairConstraintCreator::NamedAtomPairConstraintCreator() {}
NamedAtomPairConstraintCreator::~NamedAtomPairConstraintCreator() {}

core::scoring::constraints::ConstraintOP
NamedAtomPairConstraintCreator::create_constraint() const
{
	return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::NamedAtomPairConstraint( core::id::NamedAtomID(), core::id::NamedAtomID(), NULL) );
}

std::string
NamedAtomPairConstraintCreator::keyname() const {
	return "AtomPair";
}

} //namespace constraints_additional
} //namespace protocols
