// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/netcharge_energy/NetChargeConstraint.fwd.hh
/// @brief Forward declarations for a constraint for constraining sequences to have a desired net charge, analogous to a geometric constraint.
/// @details The corresponding energy term for this constraint is the NetChargeEnergy (netcharge in wts files).
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_scoring_netcharge_energy_NetChargeConstraint_fwd_hh
#define INCLUDED_core_scoring_netcharge_energy_NetChargeConstraint_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace netcharge_energy {

class NetChargeConstraint;

typedef utility::pointer::shared_ptr< NetChargeConstraint > NetChargeConstraintOP;
typedef utility::pointer::shared_ptr< NetChargeConstraint const > NetChargeConstraintCOP;

} // netcharge_energy
} // scoring
} // core

#endif

