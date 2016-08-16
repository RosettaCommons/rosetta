// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/packer_neighbors.cc
/// @brief  creates a graph that describes the possible connectivity induced by designing-in larger side chains
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/packer_neighbors.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/PackerTask.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility Headers
#include <utility/vector1.functions.hh>

#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

namespace core {
namespace pack {


/// @brief Constructs a graph where edges represent the possibility of interactions between
/// residue pairs
graph::GraphOP
create_packer_graph(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task
)
{
	utility::vector1< Distance > residue_radii = find_residue_max_radii( pose, task );

	return create_packer_graph( pose, scfxn, task, pose.total_residue(), residue_radii );
}

/// @brief Constructs a graph where edges represent the possibility of interactions between
/// residue pairs
graph::GraphOP
create_packer_graph(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	core::Size total_nodes,
	utility::vector1< Distance > const & residue_radii
)
{
	using namespace graph;

	//GraphOP g = new Graph( pose.total_residue() );
	GraphOP g( new Graph( total_nodes ) );

	if ( ! task->design_any() && ! pose.conformation().structure_moved() /* && ! core::pose::symmetry::is_symmetric( pose ) */ ) {
		//if ( false ) {
		g->copy_connectivity( pose.energies().energy_graph() );
	} else {

		/// OK -- rewriting this function to be symmetry aware
		conformation::symmetry::SymmetryInfoCOP symm_info;
		if ( dynamic_cast< conformation::symmetry::SymmetricConformation const * > ( & pose.conformation() ) ) {
			conformation::symmetry::SymmetricConformation const & symmconf(
				static_cast< conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ));
			symm_info = symmconf.Symmetry_Info();
		}

		// find radii for residues...   NOTE: flo oct 08, not anymore, doing this in above function now
		//utility::vector1< Distance > residue_radii = find_residue_max_radii( pose, task );

		// find max radius
		Distance const max_radius = utility::max( residue_radii );
		Distance const atomic_itxn_dist = scfxn.info()->max_atomic_interaction_distance();

		// create point graph and detect neighbors
		core::conformation::PointGraphOP point_graph( new core::conformation::PointGraph );
		core::conformation::residue_point_graph_from_conformation( pose.conformation(), *point_graph );
		if ( symm_info ) {
			core::conformation::find_neighbors_restricted<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( point_graph, atomic_itxn_dist + 2 * max_radius, symm_info->independent_residues() );
		} else {
			core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( point_graph, atomic_itxn_dist + 2 * max_radius );
		}

		//if ( pose.total_residue() > 468 ) {
		// Vector v316 = point_graph->get_vertex( 316 ).data().xyz();
		// Vector v468 = point_graph->get_vertex( 316 ).data().xyz();
		// std::cout << "point graph data: " << v316.x() << ", " << v316.y() << ", " << v316.z() << " and "
		//  << v468.x() << ", " << v468.y() << ", " << v468.z() << " sqrdist: " << v316.distance_squared( v468 ) << std::endl;
		//}

		// add edges
		//for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size ii = 1; ii <= total_nodes; ++ii ) {
			Size ii_asu = ii;
			if ( symm_info && symm_info->bb_follows( ii ) != 0 ) {
				ii_asu = symm_info->bb_follows(ii);
			}
			Distance const ii_itxn_rad = residue_radii[ ii_asu ] + atomic_itxn_dist;
			for ( core::conformation::PointGraph::UpperEdgeListConstIter
					iter = point_graph->get_vertex( ii ).const_upper_edge_list_begin(),
					iter_end = point_graph->get_vertex( ii ).const_upper_edge_list_end();
					iter != iter_end; ++iter ) {

				Size jj = iter->upper_vertex();
				Size jj_asu = jj;
				if ( symm_info && symm_info->bb_follows( jj ) != 0 ) {
					jj_asu = symm_info->bb_follows(jj);
				}
				Distance const jj_rad = residue_radii[ jj_asu ];

				//if ( ii == 316 && jj == 468 ) {
				// std::cout << "Packer neighbor graph creation; 316 and 468; " << iter->data().dsq() << std::endl;
				//}

				if ( jj_rad + ii_itxn_rad > 0 &&
						iter->data().dsq() < ( jj_rad + ii_itxn_rad )*( jj_rad + ii_itxn_rad ) ) {
					//std::cout << "packer_neighbors adding edge " << ii << " " << jj << std::endl;
					g->add_edge( ii, jj );
				} else {
					//std::cout << "packer_neighbors NOT adding edge " << ii << " " << jj << std::endl;
				}
			}
		}
	}
	return g;
}


/// @brief for each residue in the protein, finds the largest bounding sphere
/// over all allowable rotameric/chemical modifications possible given the input task.
///
utility::vector1< Distance >
find_residue_max_radii(
	pose::Pose const & pose,
	task::PackerTaskCOP the_task
)
{
	using namespace chemical;

	utility::vector1< Distance > residue_max_radii( pose.total_residue(), 0.0 );

	for ( Size ii = 1, ii_end = pose.total_residue(); ii <= ii_end; ++ii ) {
		Distance max_radius_for_res( 0.0 );
		/*
		if ( task->design_residue( ii ) ) {

		for ( ResidueTypeSet::AAsIter iter = residue_set.aas_defined_begin(),
		eiter = residue_set.aas_defined_end(); iter != eiter; ++iter ) {

		if ( task->allow_aa( ii, *iter )) {
		ResidueTypeCOPs const & concrete_residues( residue_set.aa_map( *iter ) );
		for ( ResidueTypeCOPs::const_iterator resiter = concrete_residues.begin(),
		eresiter = concrete_residues.end(); resiter != eresiter; ++resiter ) {
		if ( task->allow_concrete( pose.residue(ii), **resiter )) {
		if ( (*resiter)->nbr_radius() > max_radius_for_res ) {
		max_radius_for_res = (*resiter)->nbr_radius();
		}
		}
		}
		}
		}
		} else if ( task->pack_residue( ii ) ) {

		chemical::AA const resaa( pose.residue( ii ).aa());
		ResidueTypeCOPs const & concrete_residues( residue_set.aa_map( resaa ) );
		for ( ResidueTypeCOPs::const_iterator resiter = concrete_residues.begin(),
		eresiter = concrete_residues.end(); resiter != eresiter; ++resiter ) {
		if ( task->repacking_allow_concrete( pose.residue(ii), **resiter)) {
		if ( (*resiter)->nbr_radius() > max_radius_for_res ) {
		max_radius_for_res = (*resiter)->nbr_radius();
		}
		}
		}
		} */
		if ( the_task->pack_residue( ii ) ) {
			for ( task::ResidueLevelTask::ResidueTypeCOPListConstIter
					allowed_iter = the_task->residue_task( ii ).allowed_residue_types_begin(),
					allowed_end = the_task->residue_task( ii ).allowed_residue_types_end();
					allowed_iter != allowed_end; ++allowed_iter ) {
				if ( (*allowed_iter)->nbr_radius() > max_radius_for_res ) {
					max_radius_for_res = (*allowed_iter)->nbr_radius();
				}
			}
			//check whether the radius at any position needs to be increased
			Distance max_rad_change(0.0);
			for ( rotamer_set::RotSetOperationListIterator
					rotsetop_iter = the_task->residue_task( ii ).rotamer_set_operation_begin(),
					rotsetop_end = the_task->residue_task( ii ).rotamer_set_operation_end();
					rotsetop_iter != rotsetop_end; ++rotsetop_iter ) {

				core::Real radius_change = (*rotsetop_iter)->increase_packer_residue_radius( pose, the_task, ii );
				if ( radius_change > max_rad_change ) {
					max_rad_change = radius_change;
				}
			}
			if ( max_rad_change != 0.0 ) max_radius_for_res = max_radius_for_res + max_rad_change;
		} else {
			max_radius_for_res = pose.residue( ii ).nbr_radius();
		}
		residue_max_radii[ ii ] = max_radius_for_res;
	}
	return residue_max_radii;
}

/// @details pose and score function must have met before packing
/// may begin; this function will force a score evaluation if the
/// energie's scorefunction-info object does not match that of the
/// given score function.
void
pack_scorefxn_pose_handshake(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn
)
{
	// if ( true ){//pose.energies().get_scorefxn_info() != *(scfxn.info() ) ) {
	scfxn( pose );
	// }
}

} // namespace core
} // namespace pack
