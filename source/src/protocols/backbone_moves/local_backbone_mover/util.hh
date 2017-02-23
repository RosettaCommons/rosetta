// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/util.hh
/// @brief Utility functions for the local_backbone_mover.
/// @author xingjiepan (xingjiepan@gmail.com)

#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_util_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_util_hh

// Numeric headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

// Project headers
#include <protocols/backbone_moves/local_backbone_mover/types.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

/// @brief Get the coordinate of an atom given three
/// preceeding atoms and the torsion phi, angle theta
/// and distance d. Note that phi and theta should be
/// in radians.
numeric::xyzVector<Real> xyz_from_internal_coords( numeric::xyzVector<Real> atom1_xyz,
	numeric::xyzVector<Real> atom2_xyz, numeric::xyzVector<Real> atom3_xyz,
	Real phi, Real theta, Real d);


/// @brief covert an xyzVector to a vector of 3 real numbers
void xyz_to_vec1(numeric::xyzVector<Real> const& vec_xyz, vector1<Real> &vec1);


} //protocols
} //backbone_moves
} //local_backbone_mover


#endif //protocols/backbone_moves/local_backbone_mover_util_hh

