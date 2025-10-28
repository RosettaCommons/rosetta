// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/clusters/CDRCluster.hh
/// @brief Simple class to hold and access CDRcluster info at a region in the pose.  Construct in CDRClusterSet.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRCLUSTER_HH
#define INCLUDED_protocols_antibody_clusters_CDRCLUSTER_HH

#include <protocols/antibody/clusters/CDRCluster.fwd.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/antibody/AntibodyEnum.hh> // AUTO IWYU For CDRNameEnum
#include <utility/VirtualBase.hh> // AUTO IWYU For VirtualBase

#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace antibody {
namespace clusters {


/// @brief Simple class to hold and access CDRCluster info for a region of the pose.
class CDRCluster : public utility::VirtualBase {
public:
	CDRCluster(
		core::pose::Pose const & pose,
		CDRNameEnum const cdr,
		core::Size const cdr_length,
		CDRClusterEnum const cluster,
		core::Size const start,
		core::Real const distance,
		bool cis_trans_match = true);

	CDRCluster(CDRCluster const & src);

	CDRClusterOP clone() const;

	~CDRCluster() override;

	CDRNameEnum
	cdr() const { return cdr_; }

	CDRClusterEnum
	cluster() const { return cluster_; }


	/// @brief return Rosetta start that was used for construction
	core::Size
	start() const { return start_; }

	/// @brief return Rosetta end that was used for construction
	core::Size
	end() const { return end_; }


	/// @brief return PDB start.  Useful for CDR length changes in other parts of the antibody when combined with numbering scheme ala AbInfo
	core::Size
	pdb_start() const { return pdb_start_; }

	/// @brief return PDB end.  Useful for CDR length changes in other parts of the antibody when combined with numbering scheme ala AbInfo
	core::Size
	pdb_end() const { return pdb_end_; }

	/// @brief return PDB chain
	std::string
	chain() const {return chain_;}


	core::Real
	distance() const { return distance_; }

	core::Real
	length_normalized_distance() const { return normalized_distance_; }

	core::Real
	normalized_distance_in_degrees() const;

	/// @brief Does the closest cluster match at cis_trans positions?
	/// Currently, this should always be True.
	bool
	cis_trans_match() const { return cis_trans_match_; }

private:

	void
	set_pdb_numbering(core::pose::Pose const & pose, core::Size start, core::Size end);

private:

	CDRNameEnum cdr_;
	CDRClusterEnum cluster_;

	core::Real distance_;
	core::Real normalized_distance_;

	core::Size pdb_start_;
	core::Size pdb_end_;
	char pdb_start_insertion_code_;
	char pdb_end_insertion_code_;

	core::Size start_;
	core::Size end_;

	std::string chain_;

	bool cis_trans_match_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CDRCluster();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_antibody_clusters_CDRCluster )
#endif // SERIALIZATION


#endif //#ifndef INCLUDED_protocols/antibody_design_CDRCLUSTER_HH

