// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CheckForSandwichFeatures.hh
/// @brief CheckForSandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_CheckForSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_CheckForSandwichFeatures_HH

//Devel
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>
//#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFragment.hh>

using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {


	std::pair<float, float>
	cal_min_avg_dis_between_two_sheets_by_cen_res (
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size sheet_id_1,
		core::Size sheet_id_2);


	float
	cal_min_dis_between_sheets_by_all_res (
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		utility::vector1<core::Size>	all_distinct_sheet_ids);


	float
	cal_min_dis_between_two_sheets_by_all_res (
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size sheet_id_1,
		core::Size sheet_id_2);

	bool
	check_whether_sw_by_sh_id_still_alive(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id);


	std::vector<Size>
	get_all_residues_in_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sheet_id);

	std::vector<Size>
	get_cen_residues_in_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sheet_id);


	core::Real
	get_closest_distance_between_strands(
		core::pose::Pose const & pose,
		SandwichFragment strand_i,
		SandwichFragment strand_j);

	utility::vector1<SandwichFragment>
	get_full_strands_from_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);

	core::Size
	get_num_strands_in_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);


	std::string
	get_sheet_antiparallel_info(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);

	float round_to_float(
		float x);

	core::Real round_to_Real(
		core::Real x);


	core::Size round_to_Size(
		core::Real x);


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* CheckForSandwichFeatures_HH_ */
