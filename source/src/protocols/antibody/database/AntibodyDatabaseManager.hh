// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodyDatabaseManager.hh
/// @brief Handles all loading of CDR, Framework, and cluster/dmap info from external database file.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_database_AntibodyDatabaseManager_hh
#define INCLUDED_protocols_antibody_database_AntibodyDatabaseManager_hh


#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.fwd.hh>
//#include <protocols/antibody/design/GeneralAntibodyModeler.hh>
//#include <protocols/antibody/design/AntibodyDesignMover.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/AntibodyInfo.hh>

//Core Headers
#include <core/pose/Pose.hh>

//C++ Headers
#include <cppdb/frontend.h>
#include <string>
#include <map>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

namespace protocols {
namespace antibody {

//Structure for holding loaded CDRs, as well as associated information.
//Could be its own light class.
struct CDRDBPose {
	core::pose::PoseOP pose;

	clusters::CDRClusterEnum cluster;
	std::string pdb;
	core::Real distance;
	core::Real normalized_distance;
	core::Real resolution;
	//std::string germline;
	//std::string species;
	//core::Real dmap_coord;
	//species


	//Constructor
	CDRDBPose(core::pose::PoseOP pose, clusters::CDRClusterEnum cluster, std::string pdb, core::Real distance, core::Real normalized_distance, core::Real resolution) :
		pose(pose),
		cluster(cluster),
		pdb(pdb),
		distance(distance),
		normalized_distance(normalized_distance),
		resolution(resolution)
	{
	}

};

struct CDRDBSequence {
	std::string sequence;

	clusters::CDRClusterEnum cluster;
	std::string pdb;
	core::Real distance;
	core::Real normalized_distance;
	core::Real resolution;
	//std::string germline;
	//std::string species;
	//core::Real dmap_coord;
	//species


	//Constructor
	CDRDBSequence(std::string sequence, clusters::CDRClusterEnum cluster, std::string pdb, core::Real distance, core::Real normalized_distance, core::Real resolution) :
		sequence(sequence),
		cluster(cluster),
		pdb(pdb),
		distance(distance),
		normalized_distance(normalized_distance),
		resolution(resolution) {}

};

typedef std::map<CDRNameEnum, utility::vector1<CDRDBPose> > CDRDBPoseSet;
typedef std::map<CDRNameEnum, utility::vector1<CDRDBSequence> > CDRDBSequenceSet;

typedef std::map< core::chemical::AA, core::Real > AAProbabilities;




/// @brief Class responsible for loading information from the AntibodyDatabase.
class AntibodyDatabaseManager : public utility::pointer::ReferenceCount{
public:

	AntibodyDatabaseManager(AntibodyInfoCOP ab_info, bool force_north_paper_db = false);

	AntibodyDatabaseManager(AntibodyInfoCOP ab_info, std::string const database_path);

	virtual ~AntibodyDatabaseManager();

	/// @brief Set to parse outliers for both sequence and structural data.  Default False.
	void
	set_outlier_use(bool use_outliers);

	/// @brief Load CDRs from options of all cdrs..
	/// @details Will attempt to use Pose Datacache for identification of CDR clusters first.
	///
	CDRDBPoseSet
	load_cdr_poses(
		AntibodyCDRSetOptions const & options,
		core::pose::Pose const & pose,
		core::Size overhang=3);

	/// @brief Load CDRs from single cdr options.
	/// @details Will attempt to use Pose Datacache for identification of CDR clusters first.
	///
	//CDRDBPoseSet
	//load_cdr_poses(CDRDBPoseSetOptions const & options, core::pose::Pose const & pose, bool const use_light_chain_type = true, core::Size overhang=3);

	/// @brief Load probability data for CDR Design.  Returns CDRs where data could not be loaded. Cutoff indicates number of total sequences needed to use the data.
	/// @details Will attempt to use Pose Datacache for identification of CDR clusters first.
	///
	utility::vector1< bool >
	load_cdr_design_data(
		design::AntibodyCDRSeqDesignOptions const & options,
		core::pose::Pose const & pose,
		std::map< core::Size, AAProbabilities > & prob_set,
		core::Size const cutoff);

	utility::vector1< bool >
	load_cdr_design_data_for_cdrs(
		utility::vector1<bool> const & cdrs,
		core::pose::Pose const & pose,
		std::map< core::Size, AAProbabilities > & prob_set,
		core::Size const cutoff);

	/// @brief Load a set of CDR sequences from the database
	CDRDBSequenceSet
	load_cdr_sequences(
		utility::vector1< bool > const & cdrs,
		core::pose::Pose const & pose,
		bool match_on_length = false);

	//CDRSequenceSet
	//load_cdr_sequences(
	// utility::vector1< bool > cdrs,
	// core::pose::Pose const & pose,
	// AntibodyCDRSetOptions const & options,
	// bool use_light_chain_type = true,
	// bool match_on_length = false)

	utility::vector1< CDRDBSequence >
	load_cdr_sequences_for_cdr_length(
		CDRNameEnum const cdr,
		core::Size const length);

	utility::vector1< CDRDBSequence >
	load_cdr_sequences_for_cdr_cluster(
		CDRNameEnum const cdr,
		clusters::CDRClusterEnum const cluster);

	//std::map< core::Size, AAProbabilities >
	//load_framework_design_data(AntibodyInfoCOP ab_info, core::pose::Pose const & pose);

	/// @brief Use the light chain type for CDRPoses and CDRSequences if set in AntibodyInfo.
	void
	use_light_chain_type(bool use_light_chain_type);

	bool
	use_light_chain_type() const {
		return use_light_chain_type_;
	}

private:



	utility::vector1< CDRDBSequence >
	load_cdr_sequences_for_cdr(
		CDRNameEnum const cdr,
		core::Size const length,
		clusters::CDRClusterEnum const cluster,
		bool load_on_length);

	void
	start_database_session(std::string const database_path);


	/// @brief  Checks to make sure the instructions make sense before trying to create the statement for the db..
	void
	check_for_graft_instruction_inconsistencies(AntibodyCDRSetOptions const & options);



	/// @brief Checks for inconsistency in include_only and leave_out string vectors from GraftInstructions.
	template < typename T >
	bool
	has_vec_inconsistency(vector1<T> const &  include, vector1<T> const & leave_out) const;

	//Tried explicit function specialization here, but could not get it to work with clang.

	/// @brief Bind the values in the vector to the select statement.
	void
	bind_vec_constraint(utility::vector1< std::string> const & vec, cppdb::statement  & select_statement, core::Size & col) const;


	utility::vector1<std::string>
	get_cluster_string_vec(utility::vector1<clusters::CDRClusterEnum> const & clusters);

	/// @brief Does the work.  Not pretty.
	utility::vector1<core::Size>
	get_cluster_totals(bool use_outliers) const;

	/// @brief Main function called - deals with outliers vs H3_outliers
	utility::vector1<core::Size>
	get_cluster_totals() const;

	std::string
	get_ands(std::string name, core::Size n) const;

private:

	//protocols::features::ProteinSilentReportOP protein_silent_report_;
	std::string db_path_;
	utility::sql_database::sessionOP db_session_;
	AntibodyInfoCOP ab_info_;
	bool use_outliers_;
	bool use_h3_graft_outliers_; //H3 does not cluster well, so most of the structures are defined as 'outliers'
	bool use_only_H3_kinked_;
	bool use_light_chain_type_;
};

}
}

#endif //INCLUDED_protocols_antibody_design_AntibodyDatabaseManager_hh
