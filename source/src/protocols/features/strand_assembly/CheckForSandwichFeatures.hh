// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CheckForSandwichFeatures.hh
/// @brief Check various properties for SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_CheckForSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_CheckForSandwichFeatures_HH

//Devel
//#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFragment.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>


using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {

	using namespace core; 


core::Real
absolute_vec (numeric::xyzVector<core::Real> vector);


core::Real
calculate_dihedral_w_4_resnums(
	core::pose::Pose const & pose,
	core::Size res1_sheet_i,
	core::Size res2_sheet_i,
	core::Size res1_sheet_j,
	core::Size res2_sheet_j);

std::vector<core::Real>
cal_dis_angle_to_find_sheet( // calculate distance and angle to find sheet
	core::pose::Pose const & pose,
	core::Size res_i_0,
	core::Size res_i_1,
	core::Size res_i_2,
	core::Size res_j_0,
	core::Size res_j_1,
	core::Size res_j_2);

std::pair<core::Real, core::Real>
cal_min_avg_dis_between_sheets_by_cen_res (
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose & dssp_pose,
	utility::vector1<core::Size>	all_distinct_sheet_ids,
	core::Size	min_num_strands_in_sheet_);


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

core::Size
cal_num_of_sheets_that_surround_this_sheet (
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose & dssp_pose,
	utility::vector1<core::Size>	all_distinct_sheet_ids,
	core::Size sheet_id,
	core::Size	min_num_strands_in_sheet_,
	core::Real	inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_);

float
cal_shortest_dis_between_facing_aro_in_sw (
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	utility::vector1<core::Size>	all_distinct_sheet_ids,
	core::Size	min_num_strands_in_sheet_);


std::string
check_canonicalness_of_LR(
	core::Size loop_size,
	bool intra_sheet,
	std::string LR);

std::string
check_canonicalness_of_PA(
	core::Size loop_size,
	bool intra_sheet,
	std::string PA_by_preceding_E,
	std::string PA_by_following_E,
	core::Real check_canonicalness_cutoff_);

std::string
check_canonicalness_of_parallel_EE(
	core::Size loop_size,
	bool intra_sheet,
	std::string parallel_EE);

std::string
check_heading_direction ( // for example, positive,
	core::pose::Pose & dssp_pose,
	core::Size preceding_E,
	core::Size following_E,
	std::string check_heading_direction_by_);

bool
check_helix_existence(
	core::pose::Pose const & pose);

std::string
check_LR (
	core::pose::Pose & dssp_pose,
	core::Size preceding_E,
	core::Size following_E);

std::pair<std::string, std::string>
check_PA(
	core::pose::Pose & dssp_pose,
	core::Size residue_begin,
	core::Size residue_end);

bool
check_strand_too_closeness	(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	core::Real	min_inter_sheet_dis_CA_CA_);


core::Real
check_sw_by_dis(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell,
	core::Real	min_sheet_dis_,
	core::Real	max_sheet_dis_);

bool
check_whether_hairpin_connects_short_strand(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size start_res,
	core::Size next_start_res);


bool
check_whether_same_direction_strands_connect_two_sheets_or_a_loop(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size start_res,
	core::Size next_start_res);

bool
check_whether_sw_by_sh_id_still_alive(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id);

std::string
check_whether_sw_is_not_connected_with_continuous_atoms(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	core::Size sw_can_by_sh_id);

bool
check_whether_this_pdb_should_be_excluded (
	std::string tag); // I don't know how to correctly extract beta-sandwich from 1W8N for now

bool
check_whether_this_sheet_is_too_short(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_i);



bool
check_whether_sheets_are_connected_with_near_bb_atoms(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	core::Size sw_can_by_sh_id,
	core::Real	min_N_O_dis_between_two_sheets_,
	core::Real	min_N_H_O_angle_between_two_sheets_);

bool
check_whether_strand_i_is_in_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size segment_id);

std::vector<core::Size>
count_AA_w_direction(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_w_center_000,
	core::Size sw_can_by_sh_id,
	core::Size sheet_id,
	core::Size residue_begin,
	core::Size residue_end);

std::vector<core::Size>
count_AA_wo_direction(
	core::pose::Pose const & pose,
	core::Size residue_begin,
	core::Size residue_end);

std::string
determine_core_heading_surface_heading_by_distance(
	core::pose::Pose const & pose_w_center_000,
	core::Size	residue_num);

std::string
determine_heading_direction_by_vector
(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size sw_can_by_sh_id,
	core::Size sheet_id,
	core::Size residue_begin,
	core::Size residue_end,
	core::Size	ii // residue_number
);

core::Size
find_sheet	(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell, // if false, try to find a sheet in parallel way
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_,
	core::Real	min_C_O_N_angle_
	);

std::vector<Size>
get_all_residues_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sheet_id);



utility::vector1<SandwichFragment>
get_all_strands_in_sheet_i(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id);


std::vector<Size>
get_aro_residues_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	core::Size sheet_id);

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
	core::Size j_resnum_3,
	core::Real	min_sheet_dis_,
	core::Real	max_sheet_dis_);

core::Real
get_avg_dis_strands(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j);

std::vector<Size>
get_central_residues_in_other_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id,
	core::Size sheet_id);

std::vector<Size>
//get_cen_residues_in_this_sheet
get_central_residues_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sheet_id);

std::pair<int, int>
get_central_residues_in_each_of_two_edge_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size sheet_i,
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_);

utility::vector1<core::Size>
get_chain_B_resNum(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size	sw_can_by_sh_id);


core::Real
get_closest_distance_between_strands(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j);


std::pair<core::Size, core::Size>
get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size sw_can_by_sh_id,
	core::Size sheet_id,
	core::Size residue_begin,
	core::Size residue_end);

utility::vector1<SandwichFragment>
get_current_strands_in_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<core::Size>
get_distinct_sheet_id_from_sheet_table(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<core::Size>
get_distinct_sw_id_from_sw_by_components_table(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<SandwichFragment>
get_full_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<SandwichFragment>
get_full_strands_from_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id);

utility::vector1<Size>
get_list_of_residues_in_sheet_i(
	utility::vector1<SandwichFragment>	all_strands_in_sheet_i);


core::Size
get_max_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

std::pair<core::Size, core::Size>
get_next_starting_res_for_connecting_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id,
	core::Size former_ending_res);

core::Size
get_num_of_distinct_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

core::Size
get_num_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size	sheet_id);

core::Size
get_num_strands_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id);


core::Size
get_segment_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size all_strands_index);

std::string
get_sheet_antiparallel_info(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id);

core::Real
get_shortest_among_4_vals(
	core::Real arr_dis_inter_sheet[]);


core::Size
get_size_sw_by_components_PK_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id);

// used for judge_facing
utility::vector1<SandwichFragment>
get_start_end_res_num_in_the_longest_strand(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id);

std::pair<core::Size, core::Size>
get_starting_res_for_connecting_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id,
	core::Size former_res_end);

std::string
get_tag(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);


utility::vector1<core::Size>
get_vec_AA_kind (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id);

utility::vector1<Size>
get_vec_distinct_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size	sw_can_by_sh_id);

utility::vector1<Size>
get_vec_of_sw_can_by_sh_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<core::Size>
get_vector_of_strand_AA_distribution (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string heading_direction, // like core_heading, surface_heading
	std::string strand_location // like edge_strand, core_strand
	);

core::Size
identify_sheet_id_by_residue_end(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size residue_end);

// See whether this strand is an edge strand without 'sheet_antiparallel' info
std::string
is_this_strand_at_edge	(
	core::pose::Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size	sheet_id,
	core::Size	residue_begin,
	core::Size	residue_end,
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_);


std::string
is_this_strand_at_edge_by_looking_db(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size residue_begin);

int
judge_facing(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size sheet_i,
	core::Size sheet_j,
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_,
	core::Real	min_sheet_angle_by_four_term_cen_res_,
	core::Real	max_sheet_angle_by_four_term_cen_res_,
	core::Real	min_sheet_torsion_cen_res_,
	core::Real	max_sheet_torsion_cen_res_,
	core::Real	max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_);

std::string
report_heading_directions_of_all_AA_in_a_strand	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size sw_can_by_sh_id,
	core::Size sheet_id,
	core::Size residue_begin,
	core::Size residue_end);

utility::vector1<int>
retrieve_residue_num_of_rkde(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size	sw_can_by_sh_id,
	std::string	dssp_code,
	std::string	heading_direction);

float
round_to_float(
	float x);

core::Real
round_to_Real(
	core::Real x);

core::Size
round_to_Size(
	core::Real x);

std::string
see_edge_or_core_or_loop_or_short_edge(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size residue_num);

bool
see_whether_sheets_can_be_combined(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size i_sheet,
	core::Size j_sheet,
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_,
	core::Real	min_C_O_N_angle_);

std::string
see_whether_sheet_is_antiparallel(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size i_sheet,
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_,
	core::Real	min_C_O_N_angle_);

bool
see_whether_this_sw_has_SS_bond(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* CheckForSandwichFeatures_HH_ */
