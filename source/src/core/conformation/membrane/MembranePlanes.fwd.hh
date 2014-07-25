// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/MembranePlanes.fwd.hh
///
/// @brief 		Specification for Membrane Planes Residue Definitions
///	@details	When using the membrane code, users can optionally view the membrane planes
///				defined by the center/normal positions. This object will store the positions of anchoring residues
///				which will be used to draw CGO planes in PyMOl, representing the membrane planes.
///				This data should not be used in dynamic simulations, it's only purpose is for visualization.
///				Last Modified: 7/23/14
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembranePlanes_fwd_hh
#define INCLUDED_core_conformation_membrane_MembranePlanes_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {

class MembranePlanes;
typedef utility::pointer::owning_ptr< MembranePlanes > MembranePlanesOP;
typedef utility::pointer::owning_ptr< MembranePlanes const > MembranePlanesCOP;

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembranePlanes_fwd_hh
