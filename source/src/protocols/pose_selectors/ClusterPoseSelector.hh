// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/pose_selectors/ClusterPoseSeletor.hh
/// @brief  Cluster pose selector using protocols::cluster
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_pose_selectors_ClusterPoseSeletor_hh
#define INCLUDED_protocols_pose_selectors_ClusterPoseSeletor_hh

// Unit Headers
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>

// Project headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <vector>

namespace protocols {
namespace pose_selectors {

/// @brief Cluster poses by some real property as reported by the connected reporter (RMSD, ...)

class ClusterPoseSelector : public protocols::rosetta_scripts::PoseSelector {

public:
	ClusterPoseSelector();
	virtual ~ClusterPoseSelector() {}

	static std::string name() {	return "ClusterPoseSelector"; }
	virtual std::string get_name() const { return name(); }

	virtual
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	utility::vector1<bool> select_poses( utility::vector1< core::pose::PoseOP > poses ) const;
	
private:
	protocols::rosetta_scripts::PosePropertyReporterOP reporter_;
	core::Real radius_;
	core::Size structures_per_cluster_;
	core::Size max_cluster_size_;
	core::Size max_clusters_;
	core::Size max_structures_;
	core::Size initial_cluster_set_size_;
	bool remove_singletons_;

}; // ClusterPoseSelector

} // pose_selectors
} // protocols

#endif //INCLUDED_protocols_pose_selectors_ClusterPoseSeletor_hh
