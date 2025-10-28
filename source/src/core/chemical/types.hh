// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/types.hh
/// @brief  rosetta project type declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_core_chemical_types_hh
#define INCLUDED_core_chemical_types_hh

#include <string>

// Numeric headers

// Platform headers

// C++ headers


namespace core {
namespace chemical {

enum Hybridization {
	SP2_HYBRID = 1,
	SP3_HYBRID,
	RING_HYBRID,
	UNKNOWN_HYBRID,
	HYBRID_MAX = UNKNOWN_HYBRID
};

}

} // namespace core


#endif // INCLUDED_core_types_HH
