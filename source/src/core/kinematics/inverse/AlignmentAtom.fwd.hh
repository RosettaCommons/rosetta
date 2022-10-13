// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/inverse/AlignmentAtom.fwd.hh
/// @brief  Utility functions for calculating jumps by knowing desired atom positions
/// @author Jack Maguire

#ifndef INCLUDED_core_kinematics_inverse_AlignmentAtom_FWD_HH
#define INCLUDED_core_kinematics_inverse_AlignmentAtom_FWD_HH

namespace core {
namespace kinematics {
namespace inverse {

struct AlignmentAtom;

///@brief This class contains the 3 AlignmentAtoms that we want to use to define the new jump value.
struct AlignmentAtomArray;

} // namespace inverse
} // namespace kinematics
} // namespace core


#endif // INCLUDED_core_kinematics_inverse_AlignmentAtom_FWD_HH
