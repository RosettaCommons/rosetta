// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/clusters/CDRClusterSet
/// @brief Functions for CDRClusterEnums
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRCLUSTERSET_HH
#define INCLUDED_protocols_antibody_clusters_CDRCLUSTERSET_HH

#include <protocols/antibody/clusters/CDRClusterSet.fwd.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/CDRCluster.hh>
#include <protocols/antibody/clusters/CDRClusterMatcher.hh>

#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace antibody {
namespace clusters {

/// @brief Class that can determine a CDR Cluster, hold that information, and give it out when asked.
///
class CDRClusterSet : public utility::pointer::ReferenceCount {
public:

	/// @brief Constructor should only be used within AntibodyInfo, as it requires and is part of AbInfo.
	CDRClusterSet(AntibodyInfo * ab_info);

	// Undefined, commenting out to fix PyRosetta build  CDRClusterSet(CDRClusterSet const & src);

	virtual ~CDRClusterSet();

	/// @brief Identify the cluster of the CDR, using numbering information held in AntibodyInfo. Replace data if already present.
	void
	identify_and_set_cdr_cluster( core::pose::Pose const & pose, CDRNameEnum cdr);

	////////////////////////////////////////////////////////////////////////////
	// Data Removal
	//

	/// @brief Clear the held cluster data
	void
	clear();

	/// @brief Remove data of CDR using normal numbering.
	void
	clear(CDRNameEnum cdr);

	////////////////////////////////////////////////////////////////////////////
	// Data Access
	//

	bool
	empty() const;

	bool
	empty(CDRNameEnum cdr) const;

	CDRClusterCOP
	get_cluster_data(CDRNameEnum cdr) const;

	/// @brief Manually set the CDR Cluster
	void
	set_cluster_data(CDRNameEnum cdr, CDRClusterCOP cluster);

	CDRClusterEnum
	get_cluster(CDRNameEnum cdr)  const;


	/// @brief Get a new BasicCDRClusterSet with copies of CDRClusters contained in this set.
	BasicCDRClusterSetOP
	get_cacheable_cluster_data() const;

	/// @brief Set a new BasicCDRClusterSet to a pose with copies of CDRClusters contained in this set.
	void
	set_cacheable_cluster_data_to_pose(core::pose::Pose & pose ) const;


private:

	void
	calculate_cdr_cluster(CDRNameEnum cdr, core::Size start, core::Size end);

private:

	AntibodyInfo * ab_info_;
	CDRClusterMatcherOP cluster_matcher_;
	utility::vector1< CDRClusterOP > clusters_;

};

/// @brief Basic container class for CDRClusterSet, with some extra information
class BasicCDRClusterSet : public basic::datacache::CacheableData {
public:
	BasicCDRClusterSet();
	BasicCDRClusterSet(utility::vector1<CDRClusterOP> clusters);

	BasicCDRClusterSet(BasicCDRClusterSet const & src);

	~BasicCDRClusterSet();

	virtual basic::datacache::CacheableDataOP
	clone() const;

	/// @brief Set the CDRCluster
	void
	set_cluster(CDRNameEnum cdr, CDRClusterCOP cluster);

	/// @brief Get the CDRCluster
	CDRClusterCOP
	get_cluster(CDRNameEnum cdr) const;


	/// @brief Get the full set of CDRClusters
	void
	set_clusters( utility::vector1< CDRClusterOP > const clusters);

	/// @brief Get the full set of CDRClusters
	//utility::vector1< CDRClusterOP >
	//get_clusters() const;

private:

	/// @details
	/// One cluster per CDR as we assume the max we have is one L and one H from AntibodyInfo.
	/// The trick to have more than one LH will be to have a composite
	///  of AntibodyInfo objects and mutable AbChains.
	///
	utility::vector1< CDRClusterOP > clusters_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //clusters
} //antibody
} //protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_antibody_clusters_CDRClusterSet )
#endif // SERIALIZATION


#endif //#ifndef INCLUDED_protocols/antibody_design_CDRCLUSTERSET_HH
