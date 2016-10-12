// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/clusters/CDRClusterMatcher.hh
/// @brief Simple class for identifying CDR clusters of an antibody
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRCLUSTERMATCHER_HH
#define INCLUDED_protocols_antibody_clusters_CDRCLUSTERMATCHER_HH

#include <protocols/antibody/clusters/CDRClusterMatcher.fwd.hh>
#include <protocols/antibody/clusters/CDRCluster.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace clusters {

/// @brief Holds data for each cluster type
struct ClusterData{
	CDRNameEnum cdr;
	CDRClusterEnum cluster;

	core::Size length;
	std::string cis_trans_conf;

	utility::vector1< core::Real > phis;
	utility::vector1< core::Real > psis;
};


/// @brief Simple class for identifying CDR clusters of an antibody or protein structure.
/// @details Main interface is through CDRClusterSet/AntibodyInfo.  That is where AntibodyNumbering can be used to access specific CDRs and numbering scheme transformations.
class CDRClusterMatcher: public utility::pointer::ReferenceCount {
public:
	CDRClusterMatcher();

	virtual ~CDRClusterMatcher();

	/// @brief Get the cluster of an antibody CDR region, defined between start and end of the pose.
	/// Should it give out an OP or not?  It's a small class... I don't have any idea...
	CDRClusterOP
	get_cdr_cluster(core::pose::Pose const & pose, CDRNameEnum const cdr, core::Size start, const core::Size end) const;

	/// @brief Get the closest cluster of a region.  Used to detect CDR-like regions in normal proteins.
	CDRClusterOP
	get_closest_cluster(core::pose::Pose const & pose, core::Size const start, core::Size const end) const;


	/// @brief skip first grouping Cis and Trans for clusters in which a Cis/Trans designation currently does not exist.
	///  Default False
	bool
	allow_rama_mismatches() const {
		return allow_rama_mismatches_;
	}

	/// @brief Set to skip first grouping Cis and Trans for clusters in which a Cis/Trans designation currently does not exist.
	///  Default False
	void
	allow_rama_mismatches( bool const allow){
		allow_rama_mismatches_ = allow;
	}


private:

	void
	load_center_data();

	std::map< std::string, utility::vector1< core::Real > >
	get_pose_angles(core::pose::Pose const & pose, core::Size const start, core::Size const end) const;



private:

	std::string center_cluster_db_path_;
	utility::vector1< ClusterData > cluster_data_;

	/// @brief skip first grouping Cis and Trans for clusters in which a Cis/Trans designation currently does not exist.
	bool allow_rama_mismatches_;

};

}
}
}


#endif //#ifndef INCLUDED_protocols/antibody_design_CDRCLUSTERMATCHER_HH

