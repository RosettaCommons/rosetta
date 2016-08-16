// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ContextGraph.hh
/// @brief  Context graph type enumeration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_ContextGraphTypes_hh
#define INCLUDED_core_scoring_ContextGraphTypes_hh

namespace core {
namespace scoring {

enum ContextGraphType {
	ten_A_neighbor_graph = 1,
	twelve_A_neighbor_graph,
	centroid_neighbor_graph,
	num_context_graph_types = centroid_neighbor_graph // keep this guy last
};

} // scoring
} // core

#endif
