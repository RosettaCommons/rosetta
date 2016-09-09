// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/find_neighbors.fwd.hh
/// @brief  forward headers for find_neighbors.hh
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
///
/// @remarks Thanks to Will Sheffler for his ideas on refining this and extending it to atom neighbors
/// @remarks Adapting libRosetta code for generalized neighbor detection

#ifndef INCLUDED_core_conformation_find_neighbors_fwd_hh
#define INCLUDED_core_conformation_find_neighbors_fwd_hh

// Package Headers
#include <core/types.hh>


// Numeric headers

// ObjexxFCL headers
//#include <ObjexxFCL/KeyFArray1D.hh>
//#include <ObjexxFCL/KeyFArray2D.hh>

// Utility headers
//#include <utility/pointer/access_ptr.hh>

// boost headers

// C++ headers
#include <utility/assert.hh>
#include <limits>
#include <vector>

#include <utility/graph/UpperEdgeGraph.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzTriple.fwd.hh>


namespace core {
namespace conformation {

enum Strategy {
	NAIVE,
	AUTOMATIC,
	OCTREE,
	THREEDGRID,
	STRIPEHASH
};

// move the following typedef to top of file instead of
// within find_neighbors()
typedef numeric::xyzTriple< core::Size > CubeKey;
typedef numeric::xyzVector< core::Real > PointPosition;

template <class Vertex, class Edge>
void
find_neighbors_naive(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
);

template <class Vertex, class Edge>
void
find_neighbors_octree(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	Strategy strategy
);

template <class Vertex, class Edge>
void
find_neighbors_3dgrid(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
);

template <class Vertex, class Edge>
void
find_neighbors_naive_restricted(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection
);

template <class Vertex, class Edge>
void
find_neighbors_octree_restricted(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection,
	Strategy strategy
);


template <class Vertex, class Edge>
void
find_neighbors_3dgrid_restricted(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const &  residue_selection
);

template <class Vertex, class Edge>
core::Size
get_nearest_neighbor(
	utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Size node_id,
	core::Real neighbor_cutoff,
	Strategy strategy = AUTOMATIC
);

/*
// Commented out to make clang compile - duplication of default arguments and forward declaration of template functions confuses the compiler (and me!)
// Brian Weitzner and Sergey Lyskov 3/5/2011

template <class Vertex, class Edge>
void
find_neighbors(
utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
core::Real neighbor_cutoff,
Strategy strategy = AUTOMATIC
);

template <class Vertex, class Edge>
void
find_neighbors_restricted(
utility::pointer::shared_ptr<utility::graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
core::Real neighbor_cutoff,
utility::vector1< bool > const & residue_selection,
Strategy strategy = AUTOMATIC
);
*/


}
}

#endif
