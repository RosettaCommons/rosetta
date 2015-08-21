// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/CheckForSandwichFeatures.hh
/// @brief Check various properties for SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_CheckForSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_CheckForSandwichFeatures_HH

//Devel
//#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFragment.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace core;
using namespace std;

Real
absolute_vec (numeric::xyzVector<Real> vector);

Real
calculate_dihedral_w_4_resnums(
	core::pose::Pose const & pose,
	Size res1_sheet_i,
	Size res2_sheet_i,
	Size res1_sheet_j,
	Size res2_sheet_j);

vector<Real>
cal_dis_angle_to_find_sheet( // calculate distance and angle to find sheet
	core::pose::Pose const & pose,
	Size res_i_0,
	Size res_i_1,
	Size res_i_2,
	Size res_j_0,
	Size res_j_1,
	Size res_j_2);

pair<Real, Real>
cal_min_avg_dis_between_sheets_by_cen_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size min_num_strands_in_sheet_);


pair<float, float>
cal_min_avg_dis_between_two_sheets_by_cen_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sheet_id_1,
	Size sheet_id_2);


float
cal_min_dis_between_sheets_by_all_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids);


float
cal_min_dis_between_two_sheets_by_all_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sheet_id_1,
	Size sheet_id_2);

Size
cal_num_of_sheets_that_surround_this_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size sheet_id,
	Size min_num_strands_in_sheet_,
	Real inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_);

float
cal_shortest_dis_between_facing_aro_in_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size min_num_strands_in_sheet_);


string
check_canonicalness_of_LR(
	Size loop_size,
	bool intra_sheet,
	string LR);

string
check_canonicalness_of_PA(
	Size loop_size,
	bool intra_sheet,
	string PA_by_preceding_E,
	string PA_by_following_E,
	Real check_canonicalness_cutoff_);

string
check_canonicalness_of_parallel_EE(
	Size loop_size,
	bool intra_sheet,
	string parallel_EE);

string
check_heading_direction ( // for example, positive,
	core::pose::Pose & dssp_pose,
	Size preceding_E,
	Size following_E,
	string check_heading_direction_by_);

bool
check_helix_existence(
	core::pose::Pose const & pose);

string
check_LR (
	core::pose::Pose & dssp_pose,
	Size preceding_E,
	Size following_E);

pair<string, string>
check_PA(
	core::pose::Pose & dssp_pose,
	Size residue_begin,
	Size residue_end);

bool
check_strand_too_closeness (
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	Real min_inter_sheet_dis_CA_CA_);


Real
check_sw_by_dis(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell,
	Real min_sheet_dis_,
	Real max_sheet_dis_);

bool
check_whether_hairpin_connects_short_strand(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size start_res,
	Size next_start_res);


bool
check_whether_same_direction_strands_connect_two_sheets_or_a_loop(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size start_res,
	Size next_start_res);

bool
check_whether_sw_by_sh_id_still_alive(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

string
check_whether_sw_is_not_connected_with_continuous_atoms(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sw_can_by_sh_id);

bool
check_whether_this_pdb_should_be_excluded (
	string tag); // I don't know how to correctly extract beta-sandwich from 1W8N for now

bool
check_whether_this_sheet_is_too_short(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_i);


bool
check_whether_sheets_are_connected_with_near_bb_atoms(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sw_can_by_sh_id,
	Real min_N_O_dis_between_two_sheets_,
	Real min_N_H_O_angle_between_two_sheets_);

bool
check_whether_strand_i_is_in_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size segment_id);

vector<Size>
count_AA_w_direction(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_w_center_000,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end);

vector<Size>
count_AA_wo_direction(
	core::pose::Pose const & pose,
	Size residue_begin,
	Size residue_end);

string
determine_core_heading_surface_heading_by_distance(
	core::pose::Pose const & pose_w_center_000,
	Size residue_num);

string
determine_heading_direction_by_vector
(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end,
	Size ii // residue_number
);

Size
find_sheet (
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j,
	bool antiparalell, // if false, try to find a sheet in parallel way
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Real min_C_O_N_angle_
);

vector<Size>
get_all_residues_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);


utility::vector1<SandwichFragment>
get_all_strands_in_sheet_i(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);


vector<Size>
get_aro_residues_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sheet_id);

Real
get_avg_dis_CA_CA(
	core::pose::Pose const & pose,
	Size i_resnum,
	Size i_resnum_1,
	Size i_resnum_2,
	Size i_resnum_3,
	Size j_resnum,
	Size j_resnum_1,
	Size j_resnum_2,
	Size j_resnum_3,
	Real min_sheet_dis_,
	Real max_sheet_dis_);

Real
get_avg_dis_strands(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j);

vector<Size>
get_central_residues_in_other_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size sheet_id);

vector<Size>
//get_cen_residues_in_this_sheet
get_central_residues_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);

pair<int, int>
get_central_residues_in_each_of_two_edge_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sheet_i,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_);

utility::vector1<Size>
get_chain_B_resNum(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);


Real
get_closest_distance_between_strands(
	core::pose::Pose const & pose,
	SandwichFragment strand_i,
	SandwichFragment strand_j);


pair<Size, Size>
get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end);

utility::vector1<SandwichFragment>
get_current_strands_in_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<Size>
get_distinct_sheet_id_from_sheet_table(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<Size>
get_distinct_sw_id_from_sandwich_table(
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
	Size sheet_id);

utility::vector1<Size>
get_list_of_residues_in_sheet_i(
	utility::vector1<SandwichFragment> all_strands_in_sheet_i);


Size
get_max_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

pair<Size, Size>
get_next_starting_res_for_connecting_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size former_ending_res);

Size
get_num_of_distinct_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

Size
get_num_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);

Size
get_num_strands_in_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);


Size
get_segment_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size all_strands_index);

string
get_sheet_antiparallel_info(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);

Real
get_shortest_among_4_vals(
	Real arr_dis_inter_sheet[]);


Size
get_size_sandwich_PK_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

// used for judge_facing
utility::vector1<SandwichFragment>
get_start_end_res_num_in_the_longest_strand(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id);

pair<Size, Size>
get_starting_res_for_connecting_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size former_res_end);

string
get_tag(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);


utility::vector1<Size>
get_vec_AA_kind (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

utility::vector1<Size>
get_vec_distinct_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

utility::vector1<Size>
get_vec_of_sw_can_by_sh_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

utility::vector1<Size>
get_vector_of_strand_AA_distribution (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	string heading_direction, // like core_heading, surface_heading
	string strand_location // like edge_strand, core_strand
);

Size
identify_sheet_id_by_residue_end(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size residue_end);

// See whether this strand is an edge strand without 'sheet_antiparallel' info
string
is_this_strand_at_edge (
	core::pose::Pose const & pose,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id,
	Size residue_begin,
	Size residue_end,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_);


string
is_this_strand_at_edge_by_looking_db(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size residue_begin);

int
judge_facing(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sheet_i,
	Size sheet_j,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Real min_sheet_angle_by_four_term_cen_res_,
	Real max_sheet_angle_by_four_term_cen_res_,
	Real min_sheet_torsion_cen_res_,
	Real max_sheet_torsion_cen_res_,
	Real max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_);


void process_decoy(
	core::pose::Pose & dssp_pose,
	core::scoring::ScoreFunction const&
);

string
report_heading_directions_of_all_AA_in_a_strand (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end);

utility::vector1<int>
retrieve_residue_num_of_rkde(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	string dssp_code,
	string heading_direction);

float
round_to_float(
	float x);

Real
round_to_Real(
	Real x);

Size
round_to_Size(
	Real x);

string
see_edge_or_core_or_loop_or_short_edge(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size residue_num);

bool
see_whether_sheets_can_be_combined(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size i_sheet,
	Size j_sheet,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Real min_C_O_N_angle_);

string
see_whether_sheet_is_antiparallel(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size i_sheet,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Real min_C_O_N_angle_);

bool
see_whether_this_sw_has_SS_bond(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* CheckForSandwichFeatures_HH_ */
