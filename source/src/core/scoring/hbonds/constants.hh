// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/types.hh
/// @brief  core::scoring package type declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_core_scoring_hbonds_constants_hh
#define INCLUDED_core_scoring_hbonds_constants_hh


// Project Headers
#include <core/scoring/hbonds/types.hh>
#include <core/types.hh>

// Numeric Headers

// ObjexxFCL Headers


namespace core {
namespace scoring {
namespace hbonds {


// This gives the size of the evaluation related lookup tables defined
// in the HBondDatabase
Size const HB_EVAL_TYPE_COUNT = { (hbdon_MAX-1)*(hbacc_MAX-1)*(seq_sep_MAX-1)};

//car cutoffs defining what is a hbond
static core::Real const MAX_R = { 3.0 };
static core::Real const MIN_R = { 0.0 }; // AH distance
static core::Real const MIN_xH = { -1.0 }; // cos( radians( 180.0 - 0.0 ) )  // psi cutoff -- the fade functions enforce that out-of-range interactions are not scored.
static core::Real const MIN_xD = { 0.0 }; // cos( radians( 180.0 - 90.0 ) )  // theta cutoff
static core::Real const MAX_xH = { 1.0 }; // cos( radians( 180.0 - 180.0 ) )  // psi cutoff
static core::Real const MAX_xD = { 1.0 }; // cos( radians( 180.0 - 180.0 ) )  // theta cutoff

// chi term not yet implemented
//static core::Real const MAX_xC = { numeric::constants::d::pi_2 }; // chi cutoff

static core::Real const MIN_R2 = { MIN_R * MIN_R };
static core::Real const MAX_R2 = { MAX_R * MAX_R };

static core::Real const MAX_HB_ENERGY = { 0.0 }; // at and above this cutoff, not considered a hbond

} // namespace hbonds
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_types_HH
