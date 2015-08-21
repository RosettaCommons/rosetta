// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/params.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_params_HH
#define INCLUDED_protocols_magnesium_params_HH

#include <core/types.hh>

namespace protocols {
namespace magnesium {

core::Distance const MG_HOH_DISTANCE( 2.1 );
core::Distance const MG_LIGAND_DISTANCE_CUTOFF( 3.2 );
core::Distance const MG_V_DISTANCE( 1.0 );

} //magnesium
} //protocols

#endif
