// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/MembraneSpanConstraint.fwd.hh
///
/// @brief
/// @author Jonathan Weinstein


#ifndef INCLUDED_core_scoring_constraints_MembraneSpanConstraint_fwd_hh
#define INCLUDED_core_scoring_constraints_MembraneSpanConstraint_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace constraints {


class MembraneSpanConstraint; // fwd declaration
typedef utility::pointer::shared_ptr< MembraneSpanConstraint > MembraneSpanConstraintOP;
typedef utility::pointer::shared_ptr< MembraneSpanConstraint const > MembraneSpanConstraintCOP;


} // namespace constraints
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_constraints_MembraneSpanConstraint_FWD_HH
