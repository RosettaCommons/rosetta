// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file JumpSetup
/// @brief read jump-definition file   setups fold tree an chainbreak variants
/// loop code didn't work because fold-tree to complicated ( overlapping loops )
/// @details
/// @author Oliver Lange


#ifndef INCLUDED_core_scoring_dssp_util_hh
#define INCLUDED_core_scoring_dssp_util_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace dssp {

typedef utility::vector1< PointPosition > PointList;

void
get_CA_vectors(
	PointList const & ca1, // pass by reference, so no tricks:: 3x3
	PointList const & ca2, // pass by reference, so no tricks:: 3x3
	Vector & a,
	Vector & b,
	Vector & c
);

void
get_pairing_geometry(
	pose::Pose const& pose,
	Size const res1,
	Size const res2,
	Real& orientation,
	Real& pleating1,
	Real& pleating2
);

void
get_pleating(
	pose::Pose const& pose,
	Size const pos1,
	Size const pos2,
	Size &orientation,
	Size &pleating
);


} // dssp
} // scoring
} // core

#endif
