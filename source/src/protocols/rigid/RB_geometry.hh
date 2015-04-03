// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RB_moves - rigid body geometry
/// @brief functions that set up the geometry required for rigid body moves
/// @author Monica Berrondo
/// edited 8.05.08 by G.Lemmon.  I changed jump_ids from int to core::Size.
/// I also added 3 methods to get downstream and upstream centroids so we can avoid "dummy" centroids


#ifndef INCLUDED_protocols_rigid_RB_geometry_hh
#define INCLUDED_protocols_rigid_RB_geometry_hh


// Package headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>

// C++ Headers

#include <utility/vector1.hh>


namespace protocols {
namespace geometry {

numeric::xyzMatrix_double
random_reorientation_matrix(const double phi_range= 360.0, const double psi_range= 360.0);

void
centroids_by_jump(
	core::pose::Pose const & pose,
	core::Size const jump_id,
	core::Vector & upstream_ctrd, //< output
	core::Vector & downstream_ctrd //< output
);

void
centroids_by_jump(
	core::pose::Pose const & pose,
	core::Size const jump_id,
	core::Vector & upstream_ctrd, //< output
	core::Vector & downstream_ctrd, //< output
	utility::vector1< bool > ok_for_centroid_calculation
);


std::pair < core::Vector, core::Vector > centroid_pair_by_jump(
	core::pose::Pose const & pose,
	core::Size jump_id
);

core::Vector downstream_centroid_by_jump(
	core::pose::Pose const & pose,
	core::Size jump_id
);

core::Vector upstream_centroid_by_jump(
	core::pose::Pose const & pose,
	core::Size jump_id
);

void
centroids_by_jump_int(
	core::pose::Pose const & pose,
	core::Size jump_id,
	core::Vector & upstream_ctrd,
	core::Vector & downstream_ctrd
);


} // geometry
} // core

#endif
