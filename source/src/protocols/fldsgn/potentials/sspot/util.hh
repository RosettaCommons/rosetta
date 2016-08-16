// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief helper functions for SSPairPotential and HSPairPotential
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_util_HH
#define INCLUDED_protocols_fldsgn_potentials_sspot_util_HH

/// Unit headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

typedef core::Size Size;
typedef core::Real Real;
typedef core::Vector Vector;
typedef core::pose::Pose Pose;

void
spherical(
	Vector const & a2,
	Vector const & a4,
	Real & phi,
	Real & theta,
	Vector const & cen1,
	Vector const & cen2,
	Vector const & v21
);

/// @brief identifies the sequence separation along the fold tree
/// add the gap_size (default 10 ) into seqsep when there is chain break between pos1 and pos2
Size
get_foldtree_seqsep(
	Pose const & pose,
	Size pos1,
	Size pos2,
	Size gap_size=10
);


} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols

#endif
