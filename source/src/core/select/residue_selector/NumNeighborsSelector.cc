// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/NumNeighborsSelector.cc
/// @brief  The NumNeighborsSelector identifies all residues that have at least X neighbors within a distance D.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/residue_selector/NumNeighborsSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/pose/Pose.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

namespace core {
namespace select {
namespace residue_selector {


NumNeighborsSelector::NumNeighborsSelector() :
	count_water_( false ),
	threshold_( 17 ),
	distance_cutoff_( 10.0 )
{}
NumNeighborsSelector::~NumNeighborsSelector() {}

NumNeighborsSelector::NumNeighborsSelector( Size threshold, Real distance_cutoff ) :
	count_water_( false ),
	threshold_( threshold ),
	distance_cutoff_( distance_cutoff )
{}

ResidueSubset
NumNeighborsSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueSubset subset( pose.total_residue(), false );

	conformation::PointGraphOP pg( new conformation::PointGraph );
	conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );
	conformation::find_neighbors( pg, distance_cutoff_ );

	if ( count_water_ ) {
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pg->get_vertex( ii ).num_neighbors() >= threshold_ ) {
				subset[ ii ] = true;
			}
		}
	} else {
		utility::vector1< Size > non_water_neighbor_count( pose.total_residue(), 0 );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).aa() == chemical::aa_h2o ) continue;
			for ( conformation::PointGraph::VertexClass::UpperEdgeListIter
					iter = pg->get_vertex(ii).upper_edge_list_begin(),
					iter_end = pg->get_vertex(ii).upper_edge_list_end();
					iter != iter_end; ++iter ) {
				if ( pose.residue( iter->upper_vertex() ).aa() != chemical::aa_h2o ) {
					++non_water_neighbor_count[ ii ];
					++non_water_neighbor_count[ iter->upper_vertex() ];
				}
			}
			if ( non_water_neighbor_count[ ii ] >= threshold_ ) {
				subset[ ii ] = true;
			}
		}
	}
	return subset;
}

void NumNeighborsSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	count_water( tag->getOption< bool >( "count_water", false ));
	threshold( tag->getOption< Size >( "threshold", 17 ) );
	distance_cutoff_ = tag->getOption< Real >( "distance_cutoff", 10.0 );
}

std::string NumNeighborsSelector::get_name() const {
	return NumNeighborsSelector::class_name();
}

std::string NumNeighborsSelector::class_name() {
	return "NumNeighbors";
}

bool NumNeighborsSelector::count_water() const { return count_water_; }
Size NumNeighborsSelector::threshold() const { return threshold_; }
Real NumNeighborsSelector::distance_cutoff() const { return distance_cutoff_; }
void NumNeighborsSelector::count_water( bool setting ) { count_water_ = setting; }
void NumNeighborsSelector::threshold( Size setting ) { threshold_ = setting; }
void NumNeighborsSelector::distance_cutoff( Size setting ) { distance_cutoff_ = setting; }

ResidueSelectorOP
NumNeighborsSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new NumNeighborsSelector );
}

std::string
NumNeighborsSelectorCreator::keyname() const {
	return NumNeighborsSelector::class_name();
}

} //namespace residue_selector
} //namespace select
} //namespace core

