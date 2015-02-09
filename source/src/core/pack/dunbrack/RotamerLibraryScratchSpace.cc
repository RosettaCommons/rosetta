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
	
//	drotprob_dbb_( utility::vector1< Real >( n_bb, 0.0 ) ),
//	dneglnrotprob_dbb_( utility::vector1< Real >( n_bb, 0.0 ) ),
	
//	dchimean_dbb_( utility::vector1< Real4 >( n_bb, *new Real4 ) ),
//	dchisd_dbb_( utility::vector1< Real4 >( n_bb, *new Real4 ) ),
//	dchidevpen_dbb_( utility::vector1< Real >( n_bb, 0.0 ) ),
	
//	dE_dbb_( utility::vector1< Real >( n_bb, 0.0 ) ),
//	dE_dbb_dev_( utility::vector1< Real >( n_bb, 0.0 ) ),
//	dE_dbb_rot_( utility::vector1< Real >( n_bb, 0.0 ) ),
//	dE_dbb_semi_( utility::vector1< Real >( n_bb, 0.0 ) ),
	
//	dE_dbb_dev_perchi_( utility::vector1< Real4 >( n_bb, *new Real4 ) ),
	
	fa_dun_tot_( 0.0 ),
	fa_dun_rot_( 0.0 ),
	fa_dun_semi_( 0.0 ),
	fa_dun_dev_( 0.0 ),
	entropy_( 0.0 )
//	dentropy_dbb_( utility::vector1< Real >( n_bb, 0.0 ) )
	
{
/*	for ( Size i = 1; i <= n_bb; ++i ) {
		dchimean_dbb_.push_back( *new Real4 );
		dchisd_dbb_.push_back( *new Real4 );
	
		dE_dbb_dev_perchi_.push_back( *new Real4 );
	}*/
}

RotamerLibraryScratchSpace::~RotamerLibraryScratchSpace() {}

}
}
}

