// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/PointGraph.hh
/// @brief  Graph for detecting neighbors; vertices store points, edges store square distances
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Sam DeLuca
/// @author Doug Renfrew

#ifndef INCLUDED_core_conformation_PointGraph_hh
#define INCLUDED_core_conformation_PointGraph_hh

#ifdef PYROSETTA
	#include <core/graph/UpperEdgeGraph.hh>
	#include <utility/vector1.hh>
#endif

// Unit Headers
#include <core/conformation/PointGraph.fwd.hh>

// Package Headers
//because PointGraph is a typedef, not a 'real class', this header needs to include other full headers for downstream compilation.  Think of it as PointGraph inheriting from these two.

// Project Headers
#include <core/conformation/Conformation.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATON

namespace core {
namespace conformation {

void
residue_point_graph_from_conformation(
	Conformation const & conformation,
	PointGraph & pg
);

#ifdef PYROSETTA
	class PointGraphPy : public PointGraph {};

	/*template <class T1, class T2  >
	class UpperEdgeGraphPy : public graph::UpperEdgeGraph < T1, T2 > {   };

	typedef UpperEdgeGraphPy< PointGraphVertexData, PointGraphEdgeData > PointGraph_;

	class PointGraphPy_ : public PointGraph_ {}; */

	//typedef UpperEdgeGraphPy< PointGraphVertexData, PointGraphEdgeData > PointGraph_;

	//class PointGraphPy : public UpperEdgeGraphPy< PointGraphVertexData, PointGraphEdgeData > {};  // for PyRosetta
#endif

#ifdef    SERIALIZATION
template < class Archive > void save( Archive & arc, PointGraph const & pointgraph );
template < class Archive > void load( Archive & arc, PointGraph & pointgraph );
#endif // SERIALIZATION


} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_PointGraph )
#endif // SERIALIZATION

#endif
