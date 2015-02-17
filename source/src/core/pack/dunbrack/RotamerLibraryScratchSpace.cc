// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/RotamerLibraryScratchSpace.hh
/// @brief  Declaration of scratch space class for Dunbrack rotamer library
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <utility/vector1.hh>

// Package headers

namespace core {
namespace pack {
namespace dunbrack {


/// @details All the fixedsizearrays are allocated and initialized to 0
RotamerLibraryScratchSpace::RotamerLibraryScratchSpace() :
	utility::pointer::ReferenceCount(),
	rotprob_( 0.0 ),
	negln_rotprob_( 0.0 ),
	fa_dun_tot_( 0.0 ),
	fa_dun_rot_( 0.0 ),
	fa_dun_semi_( 0.0 ),
	fa_dun_dev_( 0.0 ),
	entropy_( 0.0 )
{
}

RotamerLibraryScratchSpace::~RotamerLibraryScratchSpace() {}

}
}
}

