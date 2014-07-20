// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/types.hh
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

// this used to be in PDB_Info.hh and had a space as the first character -- rhiju.
static std::string const chr_chains( "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz" );

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
