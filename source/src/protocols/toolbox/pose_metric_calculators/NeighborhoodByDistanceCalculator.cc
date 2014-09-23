// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NeighborhoodByDistanceCalculator.cc
/// @brief NeighborhoodByDistanceCalculator can determine all the neighbors of group of residues within a certain distance.
/// @author Steven Lewis


//Unit headers
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>

#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>

//Utility headers
#include <basic/options/option.hh>

//#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

static thread_local basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.NeighborhoodByDistanceCalculator" );
using basic::Error;
using basic::Warning;

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

///@details
NeighborhoodByDistanceCalculator::NeighborhoodByDistanceCalculator( std::set< core::Size > central_residues )
	: parent(), central_residues_(central_residues), dist_cutoff_ ( basic::options::option[basic::options::OptionKeys::pose_metrics::neighbor_by_distance_cutoff] ), num_neighbors_(0)
{

}

///@details
NeighborhoodByDistanceCalculator::NeighborhoodByDistanceCalculator( std::set< core::Size > central_residues, core::Real dist_cutoff )
	: parent(), central_residues_(central_residues), dist_cutoff_(dist_cutoff), num_neighbors_(0)
{
}

NeighborhoodByDistanceCalculator::NeighborhoodByDistanceCalculator( NeighborhoodByDistanceCalculator const & calculator ) : parent(), central_residues_(calculator.central_residues_), dist_cutoff_(calculator.dist_cutoff_)
{}

core::pose::metrics::PoseMetricCalculatorOP NeighborhoodByDistanceCalculator::clone() const
{ return core::pose::metrics::PoseMetricCalculatorOP( new NeighborhoodByDistanceCalculator(*this) ); }

void
NeighborhoodByDistanceCalculator::lookup(
  std::string const & key,
  basic::MetricValueBase * valptr
) const
{
	//emacs is having a hissy fit about the tab whitespace in this function (choking on the (<<< >>*>())->() combo), so don't freak out if there are spaces instead of tabs and/or it displays poorly...

	if ( key == "central_residues" ) {
		basic::check_cast( valptr, &central_residues_, "central_residues expects to return a std::set< core::Size >" );
		(static_cast<basic::MetricValue< std::set< core::Size > >* > (valptr))->set( central_residues_ );

	} else if ( key == "dist_cutoff" ) {
		basic::check_cast( valptr, &dist_cutoff_, "dist_cutoff expects to return a core::Real" );
		(static_cast<basic::MetricValue<core::Real> *>(valptr))->set( dist_cutoff_ );

  } else if ( key == "num_neighbors" ) {
	  basic::check_cast( valptr, &num_neighbors_, "num_neighbors expects to return a core::Size" );
		(static_cast<basic::MetricValue<core::Size> *>(valptr))->set( num_neighbors_ );
	} else if ( key == "num_neighbors_map" ) {
		basic::check_cast( valptr, &num_neighbors_map_, "num_neighbors_map expects to return a std::map<Size, Size>");
		(static_cast<basic::MetricValue< std::map< core::Size, core::Size > > *>(valptr))->set(num_neighbors_map_);
  } else if ( key == "neighbors" ) {
	  basic::check_cast( valptr, &neighbors_, "neighbors expects to return a std::set< core::Size >" );
		(static_cast<basic::MetricValue< std::set< core::Size > > *>(valptr))->set( neighbors_ );

  } else {
	  basic::Error() << "NeighborhoodByDistanceCalculator cannot compute metric " << key << std::endl;
		utility_exit();
 }

} //lookup

std::string
NeighborhoodByDistanceCalculator::print( std::string const & key ) const
{
	using namespace basic::options;
  if ( key == "central_residues" ) {
		std::string const spacer( option[ OptionKeys::packing::print_pymol_selection].value() ? "+" : " ");
		std::string res_string("");
		for( std::set< core::Size >::const_iterator it(central_residues_.begin()), end(central_residues_.end()); it != end; ++it)
			res_string += utility::to_string(*it) + spacer;
		return res_string;

	} else if ( key == "dist_cutoff" ) {
    return utility::to_string( dist_cutoff_ );

  } else if ( key == "num_neighbors" ) {
    return utility::to_string( num_neighbors_ );

	} else if ( key == "num_neighbors_map" ) {
		std::string map_string("");
		for(std::map<core::Size, core::Size>::const_iterator it(num_neighbors_map_.begin()), end(num_neighbors_map_.end()); it != end; ++it)
		map_string += utility::to_string(it->first) + ":" + utility::to_string(it->second) + " ";
		return map_string;

  } else if ( key == "neighbors" ) {
		using namespace basic::options; //this lets you get + or (space) as spacer
		std::string const spacer( option[ OptionKeys::packing::print_pymol_selection].value() ? "+" : " ");
		std::string nbrs_string("");
		for( std::set< core::Size >::const_iterator it(neighbors_.begin()), end(neighbors_.end()); it != end; ++it)
			nbrs_string += utility::to_string(*it) + spacer;
    return nbrs_string;

  }//else
  basic::Error() << "NeighborhoodByDistanceCalculator cannot compute metric " << key << std::endl;
  utility_exit();
  return "";
} //print

void
NeighborhoodByDistanceCalculator::recompute( core::pose::Pose const & pose )
{

	//clear old data
	neighbors_.clear();
	num_neighbors_map_.clear();
	num_neighbors_ = 0;

	core::Size const nres(pose.total_residue());

	//this is the expensive part!
	core::conformation::PointGraphOP pg( new core::conformation::PointGraph ); //create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); //create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, dist_cutoff_ ); //create edges
	runtime_assert(nres == pg->num_vertices());

	//PointGraph is the one-way graph, but this is inefficient for group v group calculations - we do not want to iterate over the entire graph each time.  Instead we want to visit just the nodes in one group and see if its edges are in the second group, so we need a two-way graph to prevent reiterating the lower half every time.
	core::graph::Graph neighborgraph(nres);
	for ( core::Size r(1); r <= nres; ++r){
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex(r).upper_edge_list_begin(),
						edge_end_iter = pg->get_vertex(r).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			neighborgraph.add_edge(r, edge_iter->upper_vertex());
		}
	}
	runtime_assert(nres == neighborgraph.num_nodes());
	runtime_assert(pg->num_edges() == neighborgraph.num_edges());

	//iterating through the graph is somewhat less expensive.
	//for each residue of interest
	for( std::set< core::Size >::const_iterator it(central_residues_.begin()), end(central_residues_.end()); it != end; ++it){
		//for all edges of that residue
		core::Size neighbor_counter = 0;
		for ( core::graph::Graph::EdgeListConstIter edge_iter = neighborgraph.get_node(*it)->const_edge_list_begin(),
						edge_end_iter = neighborgraph.get_node(*it)->const_edge_list_end();
					edge_iter != edge_end_iter; ++edge_iter ) {
			//the other end of this edge goes in neighbors_
			neighbors_.insert((*edge_iter)->get_other_ind(*it));
			neighbor_counter++;
		}//for all edges of a residue
		num_neighbors_map_[(*it)] = neighbor_counter;
		//a residue is its own neighbor, so this makes it part of the set
		neighbors_.insert(*it);
	}//for all central_residue_

	num_neighbors_ = neighbors_.size();

	return;
} //recompute

} //namespace pose_metric_calculators
} //namespace toolbox
} //namespace protocols
