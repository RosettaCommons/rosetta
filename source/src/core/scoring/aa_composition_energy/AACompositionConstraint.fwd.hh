// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/aa_composition_energy/AACompositionConstraint.fwd.hh
/// @brief Forward declarations for a constraint for constraining sequences to have a desired amino acid composition, analogous to a geometric constraint.
/// @details The corresponding energy term for this constraint is the AACompositionEnergy (aa_composition in wts files).
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_scoring_aa_composition_energy_AACompositionConstraint_fwd_hh
#define INCLUDED_core_scoring_aa_composition_energy_AACompositionConstraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace aa_composition_energy {

class AACompositionConstraint;

typedef utility::pointer::shared_ptr< AACompositionConstraint > AACompositionConstraintOP;
typedef utility::pointer::shared_ptr< AACompositionConstraint const > AACompositionConstraintCOP;

} // aa_composition_energy
} // scoring
} // core

#endif

