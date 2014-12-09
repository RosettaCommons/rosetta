// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Energies.cc
/// @brief  Symmetrical Energies class to store cached energies and track the residue
/// neighbor relationships
/// @author Ingemar Andre

// Unit Headers
#include <core/scoring/symmetry/SymmetricEnergies.hh>
#include <core/scoring/Energies.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

// Package Headers
#include <core/scoring/ContextGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/ContextGraphFactory.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunctionInfo.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Project Headers
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>

// Utility headers
#include <utility/exit.hh>

#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

namespace core {
namespace scoring {
namespace symmetry {

SymmetricEnergies::SymmetricEnergies() :
  Energies()
{}

/// copy ctor -- deep copy
SymmetricEnergies::SymmetricEnergies( Energies const & other ) :
  Energies( other )
{}

/// copy ctor -- deep copy
SymmetricEnergies::SymmetricEnergies( Energies & other ) :
  Energies( other )
{}

/// assignment operator -- deep copy
Energies const &
SymmetricEnergies::operator = ( Energies const & rhs )
{
  Energies::operator=( rhs );
  return *this;
}

bool
SymmetricEnergies::same_type_as_me( Energies const & other, bool recurse  /* = true */ ) const
{
   if ( ! dynamic_cast< SymmetricEnergies const * > ( &other ) ) {
      return false;
   }
   if ( recurse ) {
      return other.same_type_as_me( *this, false );
   } else {
      return true;
   }
}


///@details make a copy of this Energies( allocate actual memory for it )
EnergiesOP
SymmetricEnergies::clone() const
{
  return EnergiesOP( new SymmetricEnergies( *this ) );
}

SymmetricEnergies::~SymmetricEnergies() {}

void SymmetricEnergies::set_derivative_graph( MinimizationGraphOP dg )
{
	derivative_graph_ = dg;
}

MinimizationGraphOP SymmetricEnergies::derivative_graph()
{
	return derivative_graph_;
}

MinimizationGraphCOP SymmetricEnergies::derivative_graph() const
{
	return derivative_graph_;
}

/// @brief Add edges to the energy_graph and the context graphs according to domain map
///
/// @details Precondition: if the graph contains any edges, then all neighbor relationships between
/// pairs of residues that have not moved relative to each other are represented by the
/// presence or absence of edges in the energy graph.  If there are no edges in the
/// energy graph, then all pair inforamtion must be calculated
/// Precondition: if two residues have changed relative to one another according to
/// the domain map, then there are no edges between them.
void
SymmetricEnergies::update_neighbor_links(
	pose::Pose const & pose
)
{
	using namespace graph;
	using namespace scoring;

	//std::cout << "update_neighbor_links: interaction dist: " << scorefxn_info_->max_atomic_interaction_distance() <<
	//	std::endl;

  // find SymmInfo
  SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	conformation::PointGraphOP pg( point_graph() );
	if ( pg == 0 ) {
		pg = conformation::PointGraphOP( new conformation::PointGraph );
	}
	fill_point_graph( pose, pg );

	// According to the domain map, add some of the edges detected by the octree to
	// the energy graph and to the context graphs

	bool all_moved( energy_graph_no_state_check().num_edges() == 0 );

	utility::vector1< ContextGraphOP > context_graphs_present;
	for ( uint ii = 1, ii_end = context_graphs().size(); ii <= ii_end; ++ii ) {
		if ( context_graphs()[ ii ] ) context_graphs_present.push_back( context_graphs()[ ii ] );
	}

	for ( uint ii = 1, ii_end = pose.total_residue(); ii <= ii_end; ++ii ) {

		int const ii_map( domain_map_during_minimization(ii) );
		bool const ii_moved( ii_map == 0 || all_moved );

		Distance const iiradius( pose.residue_type( ii ).nbr_radius() );
		Distance const ii_intxn_radius( iiradius +
			get_scorefxn_info().max_atomic_interaction_distance() );

		for ( core::conformation::PointGraph::UpperEdgeListConstIter
				ii_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
				ii_end_iter = pg->get_vertex( ii ).upper_edge_list_end();
				ii_iter != ii_end_iter; ++ii_iter ) {
			uint const jj = ii_iter->upper_vertex();
			if ( ( domain_map_during_minimization(jj) != ii_map ) || ii_moved ) {

				Distance const jjradius( pose.residue_type( jj ).nbr_radius() );
				DistanceSquared const square_distance( ii_iter->data().dsq() );

				// How about we simply make sure the radii sum is positive instead of paying for a sqrt
				// if ( std::sqrt( square_distance ) < ( ii_intxn_radius + jj_res.nbr_radius() ) ) {
				if ( ii_intxn_radius + jjradius > 0 ) {
					if ( square_distance < (ii_intxn_radius + jjradius )*(ii_intxn_radius + jjradius )) {
//						bool symm_add;
//						if (SymmConf.Symmetry_Info().subunits() > 2 ) {
//							bool symm_add = symm_info.scoring_residue(jj);
//						} else {
//							symm_add = ( symm_info.bb_is_independent(jj) )
//													|| ( !symm_info.bb_is_independent(jj) &&
//																symm_info.bb_follows(jj) <= ii );
//						}
//						if ( symm_add ) {
							energy_graph_no_state_check().add_energy_edge( ii, jj, square_distance );
//						}
					}
					for ( uint kk = 1; kk <= context_graphs_present.size(); ++kk ) {
						context_graphs_present[ kk ]->conditionally_add_edge( ii, jj, square_distance );
					}
				}
			}
		}
	}
	/// Manually set the neighbour count for the energy_graph to be symmetrical
	for ( uint res = 1; res <= pose.total_residue(); ++res ) {
		if ( !SymmConf.Symmetry_Info()->fa_is_independent( res ) ) {
			int symm_res ( SymmConf.Symmetry_Info()->bb_follows( res ) );
			int neighbors_symm ( energy_graph_no_state_check().get_node( symm_res )->num_neighbors_counting_self() );
			energy_graph_no_state_check().get_node( res )->set_num_neighbors_counting_self_static( neighbors_symm );
		}
	}
	/// Manually set the neighbour count for the energy_graph to be symmetrical
	for ( uint res = 1; res <= pose.total_residue(); ++res ) {
		for ( uint kk = 1; kk <= context_graphs_present.size(); ++kk ) {
			if ( !SymmConf.Symmetry_Info()->fa_is_independent( res ) ) {
				int symm_res ( SymmConf.Symmetry_Info()->bb_follows( res ) );
				int neighbors_symm ( context_graphs_present[ kk ]->get_node( symm_res )->num_neighbors_counting_self() );
				context_graphs_present[ kk ]->get_node( res )->set_num_neighbors_counting_self_static( neighbors_symm );
			}
		}
	}
}


/// @brief determine distance cutoff threshold based on scorefxn_info_ and
/// then add edges to the PointGraph class
void
SymmetricEnergies::fill_point_graph( pose::Pose const & pose, conformation::PointGraphOP pg ) const {

	SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );

	Distance const max_pair_radius = pose::pose_max_nbr_radius( pose );
	Distance const energy_neighbor_cutoff = 2 * max_pair_radius + get_scorefxn_info().max_atomic_interaction_distance();

	Distance const context_cutoff = max_context_neighbor_cutoff();

	Distance const neighbor_cutoff = numeric::max( energy_neighbor_cutoff, context_cutoff );

	// Stuarts O( n log n ) octree algorithm
	core::conformation::find_neighbors_restricted<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff, symm_info->independent_residues() );
}

/// @brief Create a context graph.  If the requirement is external, someone other than a ScoreFunction
/// has declared that the context graph is needed (possibly by asking for it right now) so the graph
/// must also be made up-to-date.
void
SymmetricEnergies::require_context_graph_( scoring::ContextGraphType type, bool external ) const
{
	//utility::vector1< ContextGraphOP >  cgraphs( context_graphs() );
	utility::vector1< ContextGraphOP >& cgraphs( context_graphs() );
	utility::vector1< bool >& required_cgraphs( required_context_graphs() );
	assert( cgraphs[type] == 0 );
	required_cgraphs[type] = true;
	cgraphs[type] = ContextGraphFactory::create_context_graph( type );
	if ( cgraphs[type] == 0 ) {
		utility_exit_with_message( "Error: Null returned from ContextGraphFactory::create_context_graph( " + utility::to_string( type ) + ")" );
	}
	cgraphs[type]->set_num_nodes( size() );

	if ( max_context_neighbor_cutoff() < cgraphs[type]->neighbor_cutoff() ) {
		set_max_context_neighbor_cutoff( cgraphs[type]->neighbor_cutoff() );
	}

	core::pose::Pose & pose = *owner();

	if ( external ) {
		required_cgraphs[type] = true;

		using namespace graph;

		core::conformation::PointGraphOP point_graph( new core::conformation::PointGraph );
		fill_point_graph( pose, point_graph );
		for ( uint ii = 1, ii_end = size(); ii <= ii_end; ++ii ) {
			for ( core::conformation::PointGraph::UpperEdgeListConstIter
					ii_iter = point_graph->get_vertex( ii ).upper_edge_list_begin(),
					ii_end_iter = point_graph->get_vertex( ii ).upper_edge_list_end();
					ii_iter != ii_end_iter; ++ii_iter ) {
				uint const jj = ii_iter->upper_vertex();
				DistanceSquared const square_distance( ii_iter->data().dsq() );
				cgraphs[type]->conditionally_add_edge( ii, jj, square_distance );
			}
		}
	}
  SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation() ) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	/// Manually set the neighbour count for the context_graph to be symmetrical
	for ( uint res = 1; res <= pose.total_residue(); ++res ) {
		if ( !SymmConf.Symmetry_Info()->fa_is_independent( res ) ) {
			int symm_res ( SymmConf.Symmetry_Info()->bb_follows( res ) );
			int neighbors_symm ( cgraphs[ type ]->get_node( symm_res )->num_neighbors_counting_self() );
			cgraphs[ type ]->get_node( res )->set_num_neighbors_counting_self_static( neighbors_symm );
		}
	}


}


} // namespace symmetry
} // namespace scoring
} // namespace core
