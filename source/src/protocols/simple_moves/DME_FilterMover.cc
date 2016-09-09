// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TrialMover
/// @brief performs a move and accepts it according to Monte Carlo accept/reject criterion.
/// @author Monica Berrondo

#include <protocols/simple_moves/DME_FilterMover.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

// Random number generator

static THREAD_LOCAL basic::Tracer TR( "protocols.DME_FilterMover" );


#include <string>

#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

using namespace core;


/// @details  Setup a pointgraph for later use in dma calcs
conformation::PointGraphOP
setup_dme_point_graph( pose::Pose const & ref_pose, Real const threshold )
{
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( ref_pose.conformation(), *pg );
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, threshold );
	return pg;
}

/// @details  Calculate the dme using a pointgraph
Real
point_graph_dme( conformation::PointGraph const & pg, pose::Pose const & pose )
{
	Size total(0);
	Real dme(0.0);
	for ( Size i=1; i<= pose.size(); ++i ) {
		conformation::Residue const & i_rsd( pose.residue(i) );
		for ( auto
				i_iter     = pg.get_vertex( i ).const_upper_edge_list_begin(),
				i_end_iter = pg.get_vertex( i ).const_upper_edge_list_end();
				i_iter != i_end_iter; ++i_iter ) {
			Size const j = i_iter->upper_vertex();
			Real const reference_distance( std::sqrt( i_iter->data().dsq() ) );
			Real const pose_distance( i_rsd.nbr_atom_xyz().distance( pose.residue(j).nbr_atom_xyz() ) );
			dme += ( reference_distance - pose_distance ) * ( reference_distance - pose_distance );
			++total;
		}
	}
	return std::sqrt( dme / total );
}

/// @details  Keep trying to make a move with my_mover until the dme is less than our threshold or
/// max_tries is exceeded
/// @note  At the expense of a few additional pose copies we could save the pose with the best dme and use
/// that if we exceed max_tries

void
DME_FilterMover::apply( pose::Pose & pose )
{
	// setup a pointgraph for future dme calculations
	conformation::PointGraphOP pg( setup_dme_point_graph( pose, 10.0 ) );

	// for undoing moves with dme's that are too big
	pose::Pose const start_pose( pose );

	Size ntries( 0 );
	while ( true ) { // keep looping until we succeed
		++ntries;

		my_mover_->apply( pose );

		Real const dme( point_graph_dme( *pg, pose ) );

		TR.Trace << "apply: " << type() << " ntries= " << ntries << " dme= " << dme << std::endl;

		if ( ntries >= max_tries_ || dme < dme_threshold_ ) break;

		pose = start_pose;
	}
}

std::string
DME_FilterMover::get_name() const {
	return "DME_FilterMover";
}

}  // namespace moves
}  // namespace protocols
