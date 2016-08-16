// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueTypeConstraint.fwd.hh
///
/// @brief
/// @author Sarel Fleishman


#ifndef INCLUDED_core_scoring_constraints_ResidueTypeConstraint_fwd_hh
#define INCLUDED_core_scoring_constraints_ResidueTypeConstraint_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {


class ResidueTypeConstraint; // fwd declaration
typedef utility::pointer::shared_ptr< ResidueTypeConstraint > ResidueTypeConstraintOP;
typedef utility::pointer::shared_ptr< ResidueTypeConstraint const > ResidueTypeConstraintCOP;


} // namespace constraints
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_constraints_ResidueTypeConstraint_FWD_HH
