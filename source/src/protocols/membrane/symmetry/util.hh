// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/symmetry/util.hh
///
/// @brief      Quick calculations for symmetrizing membrane proteins
/// @details    Calculate membrane position based on a fully symmetrized spanning
///             topology. Part of an experiment to stabilize the symmetric scoring
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_symmetry_util_hh
#define INCLUDED_protocols_membrane_symmetry_util_hh

// Package Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/util.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

namespace protocols {
namespace membrane {
namespace symmetry {

/// @brief Symmetrize Spans
/// @details Create a spanning topology to reflect the full symmetric
/// complex instead of just the asymmetric unit
core::conformation::membrane::SpanningTopologyOP
symmetrize_spans(
	core::pose::Pose & pose,
	core::conformation::membrane::SpanningTopology & topology
);

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_symmetry_util_hh
