// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/PointGraph.fwd.hh
/// @brief  Graph for detecting neighbors; vertices store points, edges store square distances
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Sam DeLuca
/// @author Doug Renfrew


#ifndef INCLUDED_core_conformation_PointGraph_fwd_hh
#define INCLUDED_core_conformation_PointGraph_fwd_hh

// Package Headers
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {

typedef graph::UpperEdgeGraph< PointGraphVertexData, PointGraphEdgeData > PointGraph;

typedef utility::pointer::shared_ptr< PointGraph > PointGraphOP;

} // namespace conformation
} // namespace core

#endif
