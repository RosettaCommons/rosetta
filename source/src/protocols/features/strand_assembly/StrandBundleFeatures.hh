// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/StrandBundleFeatures.cc
/// @brief extract beta strand, strand pairs, sandwiches in pdb file, see wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StrandBundleFeatures for detail
/// @author Doo Nam Kim (based on Tim Jacobs' helix_assembly)
/// @overview
///  @ task 1: Identify all beta-strands
///   @ task 1-1: Write beta-strands into database
///  @ task 2: Identify all beta-sheets with these strands
///   @ task 2-1: Identify beta-sheets if their strands' two consecutive N-O pairs H-bond to each other
///   @ task 2-2: Write beta-sheets into database
///  @ task 3: Identify all beta-sandwiches with these sheets
///   @ task 3-1: Write beta-sandwiches into database

#ifndef INCLUDED_protocols_features_strand_assembly_StrandBundleFeatures_hh
#define INCLUDED_protocols_features_strand_assembly_StrandBundleFeatures_hh

//Unit
#include <protocols/features/strand_assembly/StrandBundleFeatures.fwd.hh>

//External
#include <boost/uuid/uuid.hpp>

//Protocols
#include <protocols/features/FeaturesReporter.hh>

//Devel
#include <protocols/features/strand_assembly/StrandFragment.hh>

//Utility
#include <utility/vector1.hh>

// for string return
#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

class StrandBundleFeatures : public protocols::features::FeaturesReporter
{

public:
	StrandBundleFeatures();
	void init_from_options();

	// XRW TEMP  virtual
	// XRW TEMP  std::string
	// XRW TEMP  type_name() const
	// XRW TEMP  {
	// XRW TEMP   return "StrandBundleFeatures";
	// XRW TEMP  }

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	utility::vector1<StrandFragment> get_full_strands(StructureID struct_id, utility::sql_database::sessionOP db_session);
	utility::vector1<StrandFragment> get_selected_strands(StructureID struct_id, utility::sql_database::sessionOP db_session);

	utility::vector1<StrandFragment> get_strand_pairs(StructureID struct_id, utility::sql_database::sessionOP db_session);

	utility::vector1<StrandFragment> get_strand_from_bss_id(StructureID struct_id, utility::sql_database::sessionOP db_session, core::Size bss_id);

	bool find_antiparallel (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);
	bool find_parallel (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);

	bool
	check_strand_too_closeness (
		core::pose::Pose const & pose,
		StrandFragment strand_i, StrandFragment strand_j);

	core::Real
	get_avg_dis_CA_CA(
		core::pose::Pose const & pose,
		core::Size i_resnum,
		core::Size i_resnum_1,
		core::Size i_resnum_2,
		core::Size i_resnum_3,
		core::Size j_resnum,
		core::Size j_resnum_1,
		core::Size j_resnum_2,
		core::Size j_resnum_3);

	core::Real
	check_sw_by_dis(
		core::pose::Pose const & pose,
		StrandFragment strand_i,
		StrandFragment strand_j,
		bool antiparalell);

	bool
	final_check_sw_by_dis (
		core::pose::Pose const & pose,
		StrandFragment temp_strand_ii_i,
		StrandFragment temp_strand_ii_j,
		StrandFragment temp_strand_jj_i,
		StrandFragment temp_strand_jj_j,
		bool antiparalell // if false, find parallel way
	);

	bool judge_sw_dis_too_close (
		core::pose::Pose const & pose,
		StrandFragment temp_strand_ii_i,
		StrandFragment temp_strand_ii_j,
		StrandFragment temp_strand_jj_i,
		StrandFragment temp_strand_jj_j);

	bool judge_sw_torsion (
		core::pose::Pose const & pose,
		StrandFragment temp_strand_ii_i,
		StrandFragment temp_strand_ii_j,
		StrandFragment temp_strand_jj_i,
		StrandFragment temp_strand_jj_j);

	core::Real
	judge_sw_inter_dis (
		core::pose::Pose const & pose,
		StrandFragment temp_strand_ii_i,
		StrandFragment temp_strand_ii_j,
		StrandFragment temp_strand_jj_i,
		StrandFragment temp_strand_jj_j);

	core::Size round (core::Real x);

	core::Size get_nearest_res_from_strand(
		core::pose::Pose const & pose,
		StrandFragment strand_to_be_searched,
		core::Size rounded_resnum);

	core::Real
	sheet_torsion (
		core::pose::Pose const & pose,
		StrandFragment strand_i,
		StrandFragment strand_j);

	bool judge_facing (
		core::pose::Pose const & pose,
		StrandFragment strand_ii_i,
		StrandFragment strand_ii_j,
		StrandFragment strand_jj_i,
		StrandFragment strand_jj_j);

	core::Real shortest_dis_sidechain (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);

	core::Real get_shortest_among_4 (
		core::Real val_shortest_dis_sidechain_1,
		core::Real val_shortest_dis_sidechain_2,
		core::Real val_shortest_dis_sidechain_3,
		core::Real val_shortest_dis_sidechain_4);

	core::Real sheet_dis_by_terminals (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	core::Size
		min_num_strands_to_deal_;

	core::Size
		max_num_strands_to_deal_;

	bool extract_native_only_;
	core::Size min_res_in_strand_;
	core::Size max_res_in_strand_;
	core::Real min_O_N_dis_;
	core::Real max_O_N_dis_;
	core::Real min_sheet_dis_;
	core::Real max_sheet_dis_;
	core::Real min_sheet_torsion_;
	core::Real max_sheet_torsion_;
	core::Real min_sheet_angle_;
	core::Real max_sheet_angle_;
	core::Real min_shortest_dis_sidechain_inter_sheet_;
};

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
