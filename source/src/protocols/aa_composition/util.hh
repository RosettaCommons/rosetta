// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/aa_composition/util.hh
/// @brief  Utility functions for aa_composition-related protocols.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_aa_composition_util_hh
#define INCLUDED_protocols_aa_composition_util_hh

// Unit Headers

// Scripter Headers

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Project Headers
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace aa_composition {

/// @brief Given a pose, run DSSP and populate the helices list.
/// @details Ignores helices shorter than min_length.
void find_helices_over_length( core::pose::Pose const &pose, utility::vector1< std::pair < core::Size, core::Size > > &helices, core::Size const min_length );

} //namespace aa_composition
} //namespace protocols

#endif //INCLUDED_protocols_aa_composition_util_hh
