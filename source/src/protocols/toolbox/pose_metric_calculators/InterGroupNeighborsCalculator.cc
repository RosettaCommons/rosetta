// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/InterGroupNeighborsCalculator.cc
/// @brief This is complicated, so pay attention.  This calculator is meant for finding interfaces between protein domains - like protein-protein interfaces but within a protein.  It's more flexible than that, though.  You define groups of residues within a protein (say, the N and C terminal domains).  You then define which pairs of groups you are interested in.  This calculator returns the union of the sets of residues at the interfaces between these domains/groups.  This calculator contains a superset of the functionality of some of the other calculators, but is less efficient in simple cases.  The pose does NOT have to have been scored.
/// @author Steven Lewis

//Unit headers
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.hh>

//
#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/graph/Graph.hh>

//Utility headers
//#include <basic/options/option.hh>
//#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>


//C++ headers
//#include <set>
//#include <utility>

static thread_local basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.InterGroupNeighborsCalculator" );

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

//typedef std::set< core::Size > one_group;
//typedef std::pair< onegroup > group_pair;
//typedef utility::vector1< group_pair > group_set;

InterGroupNeighborsCalculator::InterGroupNeighborsCalculator( group_set const & groups, core::Real dist_cutoff )
	: parent(), groups_(groups), dist_cutoff_(dist_cutoff), num_neighbors_(0)
		//not doing anything to std::set<core::Size> - should initialize empty
{}

InterGroupNeighborsCalculator::InterGroupNeighborsCalculator( InterGroupNeighborsCalculator const & calculator )
	: parent(), groups_(calculator.groups()), dist_cutoff_(calculator.dist_cutoff())
{}

InterGroupNeighborsCalculator::~InterGroupNeighborsCalculator() {}

core::pose::metrics::PoseMetricCalculatorOP InterGroupNeighborsCalculator::clone() const
{ return new InterGroupNeighborsCalculator(*this); }

void
InterGroupNeighborsCalculator::lookup(
  std::string const & key,
  basic::MetricValueBase * valptr
) const
{
	if ( key == "groups" ) {
		basic::check_cast( valptr, &groups_, "groups expects to return a utility::vector1< std::pair< std::set< core::Size >, std::set< core::Size > > >" );
		(static_cast<basic::MetricValue<group_set> *>(valptr))->set( groups_ );

	} else if ( key == "dist_cutoff" ) {
		basic::check_cast( valptr, &dist_cutoff_, "dist_cutoff expects to return a core::Real" );
		(static_cast<basic::MetricValue<core::Real> *>(valptr))->set( dist_cutoff_ );

	} else if ( key == "num_neighbors" ) {
		basic::check_cast( valptr, &num_neighbors_, "num_neighbors expects to return a core::Size" );
		(static_cast<basic::MetricValue<core::Size> *>(valptr))->set( num_neighbors_ );

	} else if ( key == "neighbors" ) {
		basic::check_cast( valptr, &neighbors_, "neighbors expects to return a std::set< core::Size >" );
		(static_cast<basic::MetricValue< std::set< core::Size > > *>(valptr))->set( neighbors_ );

	} else {
		basic::Error() << "InterGroupNeighborsCalculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

} //lookup

std::string
InterGroupNeighborsCalculator::print( std::string const & key ) const
{
	if ( key == "dist_cutoff" ) {
    return utility::to_string( dist_cutoff_ );

  } else if ( key == "num_neighbors" ) {
    return utility::to_string( num_neighbors_ );

  } else if ( key == "groups" || key == "neighbors" ) {
		//set up big return string for both sets
		using namespace basic::options; //this lets you get + or (space) as spacer
		std::string const spacer( option[ OptionKeys::packing::print_pymol_selection].value() ? "+" : " ");
		std::string nbrs_string("");

		if ( key == "groups" ) {
			for( core::Size i(1), vecsize(groups_.size()); i <= vecsize; ++i){
				nbrs_string += "{ (";
				for( one_group::const_iterator it(groups_[i].first.begin()), end(groups_[i].first.end()); it != end; ++it)
					nbrs_string += spacer + utility::to_string(*it);
				nbrs_string += ") ; (";
				for( one_group::const_iterator it(groups_[i].second.begin()), end(groups_[i].second.end()); it != end; ++it)
					nbrs_string += spacer + utility::to_string(*it);
				nbrs_string += ") }";
			}
			return nbrs_string;
		} else if ( key == "neighbors" ) {
			for( std::set< core::Size >::const_iterator it(neighbors_.begin()), end(neighbors_.end()); it != end; ++it)
				nbrs_string += utility::to_string(*it) + spacer;
			return nbrs_string;
		}//neighbors or groups
  }//else
  basic::Error() << "InterGroupNeighborsCalculator cannot compute metric " << key << std::endl;
  utility_exit();
  return "";
} //print

void
InterGroupNeighborsCalculator::recompute( core::pose::Pose const & pose )
{
	//clear old data
	neighbors_.clear();
	num_neighbors_ = 0;

	//Might be a good idea to error-check that all group residues are within the pose? - can assert < nres later?
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

	//iterating through the graph is somewhat less expensive.  We will need to iterate once per group pair (domain pair)
	//for each group/domain pair
	for( core::Size i(1), vecsize(groups_.size()); i <= vecsize; ++i){
		//for the first member of the group/domain, iterate through its residues
		for( one_group::const_iterator it(groups_[i].first.begin()), end(groups_[i].first.end()); it != end; ++it){
			//for all edges of that node
			for ( core::graph::Graph::EdgeListConstIter edge_iter = neighborgraph.get_node(*it)->const_edge_list_begin(),
							edge_end_iter = neighborgraph.get_node(*it)->const_edge_list_end();
						edge_iter != edge_end_iter; ++edge_iter ) {
				core::Size const other = (*edge_iter)->get_other_ind(*it);
				//at this point, *it and other are neighbors.  *it is in the "first" group, we need to see if other is in the second.
				if(groups_[i].second.find(other) != groups_[i].second.end()){
					// *it was in group 1 and other was in group 2 - store them!
					neighbors_.insert(*it);
					neighbors_.insert(other);
				} //if these are cross-group neighbors
			}//for all edges of a node
			//we also need to check if a residue is in both groups at once - it is its own neighbor, so this makes it part of the set
			if(groups_[i].second.find(*it) != groups_[i].second.end()) neighbors_.insert(*it);
		}// for all residues in a group
	}//for all group pairs

	num_neighbors_ = neighbors_.size();

	return;
} //recompute

} //namespace pose_metric_calculators
} //namespace toolbox
} //namespace protocols
