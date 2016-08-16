// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/md/MDConstants.hh
/// @brief  protocols/md package type declarations
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_md_MDConstants_hh
#define INCLUDED_protocols_md_MDConstants_hh

#include <core/types.hh>

namespace protocols {
namespace md {

const static core::Real Boltzmann( 0.83143435 );
const static core::Real GasConst( 1.9872065e-3 );
const static core::Real MDForceFactor( 418.4 ); // kcal/mole -> gAng^2/ps^2
const static core::Real MaxAccel( 3e4 ); // ~0.1Ang/fs^2
const static core::Real MaxVel( 10.0 ); // ~0.01Ang/fs

}
}

#endif
