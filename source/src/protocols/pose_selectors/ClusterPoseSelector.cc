// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/pose_selectors/ClusterPoseSelector.hh
/// @brief  Cluster pose selector using protocols::cluster
/// @author Luki Goldschmidt <lugo@uw.edu>


#ifndef INCLUDED_protocols_pose_selectors_ClusterPoseSelector_cc
#define INCLUDED_protocols_pose_selectors_ClusterPoseSelector_cc

// Unit Headers
#include <protocols/pose_selectors/ClusterPoseSelector.hh>
#include <protocols/pose_selectors/ClusterPoseSelectorCreator.hh>
#include <protocols/rosetta_scripts/PoseSelectorFactory.hh>
#include <protocols/rosetta_scripts/PosePropertyReporterFactory.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/cluster/cluster.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/sort_predicates.hh>
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

// C++ Headers
#include <string>

static thread_local basic::Tracer TR( "protocols.pose_selectors.ClusterPoseSelector" );

namespace protocols {
namespace pose_selectors {

////////////////////////////////////////////////////////////////////////
// ClusterPoseSelector

// Creator
protocols::rosetta_scripts::PoseSelectorOP ClusterPoseSelectorCreator::create_selector() const {
  return protocols::rosetta_scripts::PoseSelectorOP( new ClusterPoseSelector() );
}

// Selector
ClusterPoseSelector::ClusterPoseSelector() :
	reporter_(/* NULL */),
	radius_(0.0),
	structures_per_cluster_(1),
	max_cluster_size_(0),
	max_clusters_(0),
	max_structures_(0),
	initial_cluster_set_size_(400),
	remove_singletons_(false)
{
}

void ClusterPoseSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	// Cluster radius -- required (no sane default)
	radius_ = tag->getOption<core::Size>("radius");

	// Structures per cluster to select: default 1
	if(tag->hasOption("structures_per_cluster")) {
		structures_per_cluster_ = tag->getOption<core::Size>("structures_per_cluster");
	}

	// Defauls for below: no limit
	if(tag->hasOption("max_cluster_size")) {
		max_cluster_size_ = tag->getOption<core::Size>("max_cluster_size");
	}
	if(tag->hasOption("max_clusters")) {
		max_clusters_ = tag->getOption<core::Size>("max_clusters");
	}
	if(tag->hasOption("max_structures")) {
		max_structures_ = tag->getOption<core::Size>("max_structures");
	}

	if(tag->hasOption("initial_cluster_set_size")) {
		initial_cluster_set_size_ = tag->getOption<core::Size>("initial_cluster_set_size");
	}
	if(tag->hasOption("remove_singletons")) {
		remove_singletons_ = tag->getOption<bool>("remove_singletons");
	}

	// Children of tag are reporters
	BOOST_FOREACH( utility::tag::TagCOP const curr_tag, tag->getTags() ) {
		protocols::rosetta_scripts::PosePropertyReporterOP new_reporter(
			protocols::rosetta_scripts::PosePropertyReporterFactory::get_instance()->
				newPosePropertyReporter( curr_tag, data, filters, movers, pose )
		);
		runtime_assert( new_reporter != 0 );
		reporter_ = new_reporter;
		TR << "Defined pose property reporter of type " << curr_tag->getName() << std::endl;
		// Only first reporter used -- add warning when multiple defined?
		break;
	}
}

utility::vector1<bool> ClusterPoseSelector::select_poses(
	utility::vector1< core::pose::PoseOP > poses
)
{
	using namespace protocols::cluster;

	TR << "Applying selector " << get_name() << std::endl;

	ClusterPhilStyleOP clustering( new ClusterPhilStyle_PoseReporter(reporter_) );

	clustering->set_cluster_radius( radius_ );
	// clustering->set_population_weight( ?? );

	// Cluster initial set
	core::Size n_poses = 0;
	utility::vector1< core::pose::PoseOP >::iterator pose_it = poses.begin();

	while( pose_it != poses.end() && ++n_poses <= initial_cluster_set_size_ ) {
		core::pose::PoseOP p( *pose_it );
		clustering->apply( *p );
		++pose_it;
	}

	clustering->do_clustering( max_clusters_ > 0 ? max_clusters_ : poses.size() );
	clustering->do_redistribution();

	if(pose_it != poses.end()) {
		// Add remaining structures
		AssignToClustersMover mover_add_structures( clustering );
		while( pose_it != poses.end() ) {
			core::pose::PoseOP p( *pose_it );
			mover_add_structures.apply( *p );
			++pose_it;
		}
		clustering->do_redistribution();
	}

	// Post processing
	// clustering->sort_each_group_by_energy();
	// clustering->sort_groups_by_energy();

	if(max_cluster_size_ > 0) {
 		clustering->limit_groupsize( (int) max_cluster_size_ );
	}

	if(remove_singletons_) {
		clustering->remove_singletons();
	}

	// Report
	clustering->print_summary();
	// clustering->print_cluster_assignment();

	// Select poses from clusters
	utility::vector1<bool> selected_poses;
  selected_poses.resize(poses.size(), false);

	core::Size ncluster = 0, nstructure = 0;
	std::vector < Cluster > const & clusters = clustering->get_cluster_list();

	BOOST_FOREACH( Cluster cluster, clusters ) {
		TR.Debug << "Cluster " << ncluster << std::endl;
		for( core::Size i = 0; i < cluster.size() && i < structures_per_cluster_; ++i ) {
			TR.Debug << "  Selecting pose " << i << " => " << cluster[i]+1 << std::endl;
			selected_poses[ cluster[i]+1 ] = true; // +1 due to vector1 -- why can't we stick to 0-based indices?!
			if( ++nstructure >= max_structures_ && max_structures_ > 0 ) {
				break;
			}
		}

		++ncluster;
		if(
			(ncluster >= max_clusters_ && max_clusters_ > 0) || 
			(nstructure >= max_structures_ && max_structures_ > 0) 
		) {
			break;
		}
	}
	
	TR << "Selected " << nstructure << " poses from " << ncluster << " clusters" << std::endl;

	return selected_poses;
}

////////////////////////////////////////////////////////////////////////

} // pose_selectors
} // protocols

#endif //INCLUDED_protocols_pose_selectors_ClusterPoseSelector_cc
