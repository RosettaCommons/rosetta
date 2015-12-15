// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/PointGraph.cc
/// @brief  Graph for detecting neighbors; vertices store points, edges store square distances
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Sam DeLuca
/// @author Doug Renfrew


// Unit Headers
#include <core/conformation/PointGraph.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>


#ifdef    SERIALIZATION
// Package serialization headers
#include <core/graph/UpperEdgeGraph.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif

namespace core {
namespace conformation {

void
residue_point_graph_from_conformation(
	Conformation const & conformation,
	PointGraph & pg
)
{
	pg.set_num_vertices( conformation.size() );
	for ( platform::Size ii = 1, ii_end = conformation.size(); ii <= ii_end; ++ii ) {
		pg.get_vertex( ii ).data().xyz() = conformation.residue(ii).xyz( conformation.residue( ii ).nbr_atom() );
	}
}

#ifdef    SERIALIZATION
template < class Archive > void save( Archive & arc, PointGraph const & pointgraph )
{
	core::graph::save_to_archive( arc, pointgraph );
}

template < class Archive > void load( Archive & arc, PointGraph & pointgraph )
{
	core::graph::load_from_archive( arc, pointgraph );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( PointGraph );
#endif // SERIALIZATION

} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_REGISTER_TYPE( core::conformation::PointGraph )
CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_PointGraph )
#endif // SERIALIZATION
