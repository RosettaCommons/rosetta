// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NeighborsByDistanceCalculator.cc
/// @brief NeighborsByDistanceCalculator can determine all the neighbors of a residue within a certain distance.  The pose does not have to have been scored (have a full Energies object).  It uses the PointGraph tools to find neighbors.  There is probably a much more sophisticated way to do this with existing Graph tools but I don't know what it is.
/// @author Steven Lewis


//Unit headers
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>


#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/PointGraphData.hh>

#include <core/graph/UpperEdgeGraph.hh>

//Utility headers
//#include <basic/options/option.hh>
//#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//C++ headers
//#include <set>

static thread_local basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.NeighborsByDistanceCalculator" );

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

NeighborsByDistanceCalculator::NeighborsByDistanceCalculator( core::Size central_residue, core::Real dist_cutoff )
	: parent(), central_residue_(central_residue), dist_cutoff_(dist_cutoff), num_neighbors_(0)
		//not doing anything to std::set<core::Size> - should initialize empty
{}

NeighborsByDistanceCalculator::NeighborsByDistanceCalculator( NeighborsByDistanceCalculator const & calculator )
	: parent(), central_residue_(calculator.central_residue()), dist_cutoff_(calculator.dist_cutoff())
{}

core::pose::metrics::PoseMetricCalculatorOP NeighborsByDistanceCalculator::clone() const
{ return core::pose::metrics::PoseMetricCalculatorOP( new NeighborsByDistanceCalculator(*this) ); }

void
NeighborsByDistanceCalculator::lookup(
  std::string const & key,
  basic::MetricValueBase * valptr
) const
{
	if ( key == "central_residue" ) {
		basic::check_cast( valptr, &central_residue_, "central_residue expects to return a core::Size" );
		(static_cast<basic::MetricValue<core::Size> *>(valptr))->set( central_residue_ );

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
		basic::Error() << "NeighborsByDistanceCalculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

} //lookup

std::string
NeighborsByDistanceCalculator::print( std::string const & key ) const
{
  if ( key == "central_residue" ) {
    return utility::to_string( central_residue_ );

  } else if ( key == "dist_cutoff" ) {
    return utility::to_string( dist_cutoff_ );

  } else if ( key == "num_neighbors" ) {
    return utility::to_string( num_neighbors_ );

  } else if ( key == "neighbors" ) {
		using namespace basic::options; //this lets you get + or (space) as spacer
		std::string const spacer( option[ OptionKeys::packing::print_pymol_selection].value() ? "+" : " ");
		std::string nbrs_string("");
		for( std::set< core::Size >::const_iterator it(neighbors_.begin()), end(neighbors_.end()); it != end; ++it)
			nbrs_string += utility::to_string(*it) + spacer;
    return nbrs_string;

  }//else
  basic::Error() << "NeighborsByDistanceCalculator cannot compute metric " << key << std::endl;
  utility_exit();
  return "";
} //print

void
NeighborsByDistanceCalculator::recompute( core::pose::Pose const & pose )
{
	//clear old data
	neighbors_.clear();
	num_neighbors_ = 0;

	//if the central residue was never set, or was set outside the pose, this will cause problems!
	if( (central_residue_ < 1) || (central_residue_ > pose.total_residue())) {
		TR.Error << "central residue " << central_residue_ << " outside of pose; NBDC reporting empty set!" << std::endl;
		return;
	}

///This is not necessarily the best implementation of this - this code's real utility is that it patches cleanly into the Calculator and TaskOperation hierarchies.  If you have a better/faster implementation, replace this and feel no guilt.

	core::conformation::PointGraphOP pg( new core::conformation::PointGraph ); //create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); //create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, dist_cutoff_ ); //create edges

	for ( core::Size r(1); r <= central_residue_; ++r){
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex(r).upper_edge_list_begin(),
						edge_end_iter = pg->get_vertex(r).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			core::Size const other = edge_iter->upper_vertex();

			//we know the two nodes connected to this edge.  We want to remember the one that is NOT central_residue_, if one of them IS central_residue_.  if neither is central_residue_, do nothing.
			if (other == central_residue_) {
				neighbors_.insert(r);
			}
			else if ( r == central_residue_ ){
				neighbors_.insert(other);
			}
		}
	}

	//now account for the residue itself
	neighbors_.insert(central_residue_);
	num_neighbors_ = neighbors_.size();

	return;
} //recompute

} //namespace pose_metric_calculators
} //namespace toolbox
} //namespace protocols
