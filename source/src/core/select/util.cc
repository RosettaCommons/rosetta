// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util.cc
/// @brief Utilities to help in selecting residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Mike Tyka
/// @author Chu Wang
/// @author Daniel J. Mandell


#include <core/select/util.hh>

#include <basic/Tracer.hh>

#include <core/scoring/TenANeighborGraph.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>

#include <numeric/xyzVector.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.util" );


namespace core {
namespace select {

utility::vector1< Size >
get_residues_from_subset( utility::vector1< bool > subset, bool select){
	utility::vector1< Size > residues;
	for ( core::Size i = 1; i <= subset.size(); ++i ) {
		if ( subset[i] == select ) {
			residues.push_back( i );
		}
	}
	return residues;
}



utility::vector1< bool >
get_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & residue_positions,
	core::Real neighbor_dis
)
{
	utility::vector1< bool > selection_and_neighbors = residue_positions;

	if ( neighbor_dis <= 10.0 ) {
		fill_neighbor_residues(pose, selection_and_neighbors);

		//Make sure to turn off subset residues!
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( residue_positions[i] ) {
				selection_and_neighbors[i] = false;
			}
		}

		return selection_and_neighbors;

	} else {
		utility_exit_with_message("get_neighbor_residues only currently works for neighbor distances of 10A or less!  Please use NeighborhoodResidueSelector instead");
	}


}


void
fill_neighbor_residues(
	core::pose::Pose const & pose,
	utility::vector1< bool > & residue_positions,
	core::Real neighbor_dis
)
{
	utility::vector1< bool > selection = residue_positions;

	fill_tenA_neighbor_residues(pose, residue_positions);
	filter_neighbors_by_distance(pose, selection, residue_positions, neighbor_dis);
}


utility::vector1< bool >
trim_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & selection,
	utility::vector1<bool> const & all_neighbors,
	core::Real & dist_cutoff
)
{
	utility::vector1< bool > trimmed_neighbors = all_neighbors;
	filter_neighbors_by_distance(pose, selection, trimmed_neighbors, dist_cutoff);
	return trimmed_neighbors;

}


void
filter_neighbors_by_distance(
	core::pose::Pose const & pose,
	utility::vector1<bool> const & selection,
	utility::vector1<bool> & selection_and_neighbors,
	core::Real & dist_cutoff
)
{

	for ( Size i = 1; i <= selection_and_neighbors.size(); ++i ) {
		if ( selection_and_neighbors[ i ] == false ) continue;

		selection_and_neighbors[ i ] = false; //Get ready to change this.
		for ( Size x = 1; x <= selection.size(); ++x ) {
			if ( ! selection[x] ) continue;

			// Get the atom vectors for loop and scaffold CB, or CA if GLY
			numeric::xyzVector< Real > neighbor_vec;
			numeric::xyzVector< Real > select_vec;
			neighbor_vec = pose.residue( i ).xyz( pose.residue( i ).nbr_atom() );
			select_vec = pose.residue( x ).xyz( pose.residue( x ).nbr_atom() );
			// only keep as neighbor if dist within cutoff
			Real dist = neighbor_vec.distance( select_vec );
			if ( dist <= dist_cutoff ) {
				selection_and_neighbors[ i ] = true;
			}
		}
	}


}

utility::vector1< bool >
get_tenA_neighbor_residues(
	pose::Pose const & pose,
	utility::vector1<bool> const & residue_positions
)
{
	utility::vector1< bool > positions_and_neighbors = residue_positions;
	fill_tenA_neighbor_residues(pose, positions_and_neighbors);
	return positions_and_neighbors;
}


void fill_tenA_neighbor_residues(
	pose::Pose const & pose,
	utility::vector1<bool> & residue_positions
)
{
	//make a local copy first because we will change content in residue_positions
	utility::vector1<bool> local_residue_positions = residue_positions;
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	for ( Size i=1; i <= local_residue_positions.size(); ++i ) {
		if ( ! local_residue_positions[i] ) continue;
		utility::graph::Node const * current_node( tenA_neighbor_graph.get_node(i)); // find neighbors for this node
		for ( utility::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
				it != current_node->const_edge_list_end(); ++it ) {
			Size pos = (*it)->get_other_ind(i);
			residue_positions[ pos ] = true;
		}
	}
}

} //core
} //select


