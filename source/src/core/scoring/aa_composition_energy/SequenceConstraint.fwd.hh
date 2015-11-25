// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/aa_composition_energy/SequenceConstraint.fwd.hh
/// @brief Forward declarations for base class for constraining sequences, analogous to a geometric constraint.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_scoring_aa_composition_energy_SequenceConstraint_fwd_hh
#define INCLUDED_core_scoring_aa_composition_energy_SequenceConstraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace aa_composition_energy {

class SequenceConstraint;

typedef utility::pointer::shared_ptr< SequenceConstraint > SequenceConstraintOP;
typedef utility::pointer::shared_ptr< SequenceConstraint const > SequenceConstraintCOP;

} // aa_composition_energy
} // scoring
} // core

#endif

