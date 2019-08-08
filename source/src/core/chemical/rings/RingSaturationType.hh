// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/rings/RingSaturationType.hh
/// @brief   Enumeration definition for RingSaturationType.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_rings_RingSaturationType_HH
#define INCLUDED_core_chemical_rings_RingSaturationType_HH

namespace core {
namespace chemical {
namespace rings {

/// @brief:  An enumeration of saturation types for ring systems.
/// @details The value labels if a ring is aliphatic, aromatic, or some other saturation state.
/// @note    Currently, the only options are ALIPHATIC and AROMATIC.  In the future, other states will be added.
enum RingSaturationType {
	ALIPHATIC = 0,
	AROMATIC
};

}  // namespace rings
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_rings_RingSaturationType_HH
