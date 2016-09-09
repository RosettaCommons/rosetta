// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/Graph.hh
/// @brief  generic graph class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_graph_graph_util_hh
#define INCLUDED_utility_graph_graph_util_hh

// Package Headers

// Utility Headers

#include <utility/graph/Graph.fwd.hh>
#include <utility/vector1.hh>


// C++ headers
//#include <utility> // for std::pair?

namespace utility {
namespace graph {

/// @brief returns a vector1 of connected component descriptions:
/// each entry holds the connected-component size
/// and a representative vertex from that connected component.
/// O( V+E ).
utility::vector1< std::pair< platform::Size, platform::Size > >
find_connected_components( Graph const & g );

void
delete_all_intragroup_edges(
	Graph & g,
	utility::vector1< platform::Size > const & node_groups
);

}
}

#endif
