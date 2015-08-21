// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file .../CDRDihedralConstraintMover.hh
/// @brief Add CDR cluster-based Dihedral CircularHarmonic constraints or general CircHarmonic csts to a CDR.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_constraints_CDRDihedralConstraintMover_hh
#define INCLUDED_protocols_antibody_constraints_CDRDihedralConstraintMover_hh

#include <protocols/antibody/constraints/CDRDihedralConstraintMover.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>


#include <core/pose/Pose.hh>

#include <basic/datacache/DataMap.fwd.hh>


// Forward
namespace protocols {
namespace antibody {
namespace constraints {

///@brief Add Cluster or General Dihedral CircularHarmonic constraints to a CDR.
///  Cluster constraints currently only work for AHO renumbered CDRs.
///   (This will be rafactored to create constraints on-the-fly from cluster Mean/SD instead of from cst files.)
///
class CDRDihedralConstraintMover : public protocols::moves::Mover {
public:

	CDRDihedralConstraintMover();

	CDRDihedralConstraintMover(AntibodyInfoCOP ab_info);

	CDRDihedralConstraintMover(AntibodyInfoCOP ab_info, CDRNameEnum cdr);

	CDRDihedralConstraintMover(CDRDihedralConstraintMover const & src);

	virtual~CDRDihedralConstraintMover();

	virtual void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	);



	///@brief Attempt to add cluster-based dihedral constraints.  If this is set to false will simply add General Dihedral constraints instead.
	/// Default True
	void
	set_use_cluster_csts(bool use_cluster_csts);

	///@brief If we are set to use cluster csts and:
	///  1) the cluster is NA,
	///  2) there is sparse data for the cluster, or
	///  3) The CDR is H3 and we have chosen not to use cluster data for H3 (default),
	///    THEN add general dihedral constraints instead.
	// Default True
	void
	set_use_general_csts_on_cluster_failure(bool use_general_csts_on_failure);


	void
	set_cdr(CDRNameEnum cdr);

	virtual void
	apply(core::pose::Pose & pose);

public:

	///@brief Do not use AntibodyInfo to for cluster - use this cluster instead
	void
	set_force_cluster(clusters::CDRClusterEnum cluster);

	///@brief Remove any forced cluster settings.
	void
	set_remove_any_set_forced_cluster();

	void
	set_cluster_csts_data_cutoff(core::Size cutoff);

	///@brief Use constraints which have the means as the actual cluster means.
	///  Setting this to false will use constraints that have the cst means set as cluster center data.
	void
	set_cluster_csts_use_mean_cst_data(bool use_mean_cst_data);

	void
	set_cluster_csts_use_outlier_data(bool use_outlier_data);

	///@brief Set to use H3 cluster data for constraints if we are doing cluster-based constraints.
	/// Default False - H3 does not cluster well.  If use_general_data_on_failure is false, we will skip H3.
	void
	set_use_cluster_for_H3(bool use_cluster_for_H3);

public:
	void
	set_general_phi_sd(core::Real phi_sd);

	void
	set_general_psi_sd(core::Real psi_sd);

	///@brief By default, if cluster information is present in the datacache, we attempt to use that first.
	/// Override this behavior by setting this option to true.
	void
	set_ignore_pose_datacache(bool ignore_pose_datacache);
public:

	std::string
	get_name() const;



	protocols::moves::MoverOP
	clone() const;

	//CDRDihedralConstraintMover & operator=(CDRDihedralConstraintMover const & src);

	virtual moves::MoverOP fresh_instance() const;

private:

	void
	set_defaults();

	void
	read_command_line_options();

private:

	/// @brief Adds a harmonic constraint to a Pose CDR based on cluster type
	/// @details Currently requires North_AHO numbering.
	bool
	add_harmonic_cluster_constraint(core::pose::Pose & pose, clusters::CDRClusterEnum const cluster);

	core::Size
	get_number_of_struct_used_for_csts(clusters::CDRClusterEnum const cluster);

	/// @brief Gets the cluster constraint name.  Returns NA if not found.
	std::string
	get_harmonic_cluster_constraint_filename(clusters::CDRClusterEnum const cluster);

	std::string
	get_harmonic_cluster_constraint_db_directory();


private:

	AntibodyInfoCOP ab_info_;
	CDRNameEnum cdr_;
	std::string db_base_path_;

	bool cdr_is_set_;

	clusters::CDRClusterEnum forced_cluster_;
	bool force_cluster_;

	bool use_cluster_csts_;
	bool use_outliers_;
	bool use_mean_cst_data_;
	bool use_general_csts_on_failure_;
	bool use_cluster_for_H3_;
	bool ignore_pose_datacache_;
	core::Size cluster_data_cutoff_;
	core::Real general_phi_sd_;
	core::Real general_psi_sd_;
};

} //constraints
} //antibody
} //protocols




#endif //INCLUDED_protocols_antibody_constraints_CDRDihedralConstraintMover







