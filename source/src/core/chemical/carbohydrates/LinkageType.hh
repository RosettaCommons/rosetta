// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/LinkageType.hh
/// @brief   Enumerator definition for LinkageType.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_carbohydrates_LinkageType_HH
#define INCLUDED_core_chemical_carbohydrates_LinkageType_HH

#include <core/id/types.hh>

namespace core {
namespace chemical {
namespace carbohydrates {

/// @brief    Labels for the linkage type of the carbohydrate phi or psi angle.
/// @details  The CHI energy functions depend on the type of linkage.\n
/// Parameters for the Gaussian functions that compose the CHI energy functions are stored in a vector indexed by these
/// labels.
enum LinkageType {
	FIRST_LINK_TYPE = 1,
	ALPHA_LINKS = 1,     // used to describe alpha linkages (phi)
	BETA_LINKS,          // used to describe beta linkages (phi)
	_2AX_3EQ_4AX_LINKS,  // used to describe ->2-axial, ->3-equatorial, or ->4-axial linkages (psi)
	_2EQ_3AX_4EQ_LINKS,  // used to describe ->2-equatorial, ->3-axial, or ->4-equatorial linkages (psi)
	N_LINK_TYPES = _2EQ_3AX_4EQ_LINKS,

	LINKAGE_NA // A Null for searching.
};


//Alpha/Beta go from 180 to 180 and are for phi.
//Other two go from 0 to 360 and are for psi.


}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_LinkageType_HH
