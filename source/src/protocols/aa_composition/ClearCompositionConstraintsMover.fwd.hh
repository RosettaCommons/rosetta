// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/ClearCompositionConstraintsMover.fwd.hh
/// @brief Forward declarations for the ClearCompositionConstraintsMover.
/// @details Removes all AACompositionConstraints from a pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_aa_composition_ClearCompositionConstraintsMover_fwd_hh
#define INCLUDED_protocols_aa_composition_ClearCompositionConstraintsMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace aa_composition {

class ClearCompositionConstraintsMover;
typedef utility::pointer::shared_ptr< ClearCompositionConstraintsMover > ClearCompositionConstraintsMoverOP;
typedef utility::pointer::shared_ptr< ClearCompositionConstraintsMover const > ClearCompositionConstraintsMoverCOP;

} // aa_composition
} // protocols

#endif
