// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.hh
/// @brief Write to a DB after SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)

#ifndef INCLUDED_protocols_features_strand_assembly_WriteToDBFromSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_WriteToDBFromSandwichFeatures_HH

//Devel
#include <protocols/features/strand_assembly/CheckForSandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>


using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace core; 

bool
change_sheet_id_if_possible(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Real	min_CA_CA_dis_,
	core::Real	max_CA_CA_dis_,
	core::Real	min_C_O_N_angle_);

core::Size
delete_this_struct_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);


core::Size
delete_this_sw_can_by_sh_id_from_sw_by_comp(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id);

utility::vector1<SandwichFragment>
prepare_WriteToDB_sandwich(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

void
WriteToDB_AA_to_terminal_loops (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose & dssp_pose,
	core::Size	sandwich_PK_id_counter,
	core::Size	sw_can_by_sh_id,
	std::string tag,
	bool starting_loop,
	core::Size residue_begin,
	core::Size residue_end);



core::Size
WriteToDB_ending_loop(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose & dssp_pose,
	core::Size	sandwich_PK_id_counter,
	core::Size	sw_can_by_sh_id,
	std::string	tag,
	core::Size	max_starting_loop_size_);

core::Size
WriteToDB_long_strand_id_in_each_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id);

core::Size
WriteToDB_avg_b_factor_CB_at_each_component (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	core::Size	sw_can_by_sh_id);

core::Size
WriteToDB_dihedral_angle_between_core_strands_across_facing_sheets (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	core::Size	sw_can_by_sh_id);

core::Size
WriteToDB_dssp_ratio_in_sw(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose & dssp_pose,
	core::Size sw_can_by_sh_id);

core::Size
WriteToDB_hydrophobic_ratio_net_charge (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id);


core::Size
WriteToDB_min_avg_dis_between_sheets_by_cen_res (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id,
	core::pose::Pose & dssp_pose,
	utility::vector1<core::Size>	all_distinct_sheet_ids,
	core::Size	min_num_strands_in_sheet_);


core::Size
WriteToDB_min_dis_between_sheets_by_all_res (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id,
	core::pose::Pose & dssp_pose,
	utility::vector1<core::Size>	all_distinct_sheet_ids);

core::Size
WriteToDB_number_of_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	core::Size max_num_sw_per_pdb_,
	core::Real min_CA_CA_dis_,
	core::Real max_CA_CA_dis_);

void
WriteToDB_number_of_core_heading_W_in_sw(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id);

core::Size
WriteToDB_number_of_core_heading_LWY_in_core_strands_in_sw(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id);

void
WriteToDB_number_of_core_heading_FWY_in_sw(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id);




core::Size
WriteToDB_number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id,
	core::Size current_bs_id,
	core::Size closest_bs_id);

core::Size
WriteToDB_number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_id,
	core::Size current_bs_id,
	core::Size closest_bs_id);

core::Size
WriteToDB_number_of_edge_strands_in_each_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id);

core::Size
WriteToDB_number_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id,
	core::Size num_of_sheets_that_surround_this_sheet);

core::Size
WriteToDB_number_strands_in_each_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id);

core::Size
WriteToDB_prolines_that_seem_to_prevent_aggregation(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id);

core::Size
WriteToDB_ratio_of_core_heading_FWY_in_sw(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id,
	core::pose::Pose const & pose);

core::Size
WriteToDB_rkde(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	rkde_PK_id_counter,
	std::string tag,
	core::Size residue_number,
	std::string	residue_type);

core::Size
WriteToDB_rkde_in_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	rkde_in_strands_PK_id_counter,
	std::string tag,
	core::Size sw_can_by_sh_id,
	core::Size residue_number,
	std::string	residue_type,
	std::string	heading_direction);

core::Size
WriteToDB_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_PK_id_counter,
	core::Size sheet_id,
	core::Size segment_id);

core::Size
WriteToDB_sheet_antiparallel(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sheet_id,
	std::string antiparallel);


core::Size
WriteToDB_sheet_connectivity(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::Size sandwich_PK_id_counter,
	std::string tag,
	core::Size sw_can_by_sh_id,
	std::string loop_kind,
	core::Size intra_sheet_con_id,
	core::Size inter_sheet_con_id,
	std::string LR,
	std::string cano_LR,
	std::string	PA_by_preceding_E,
	std::string	PA_by_following_E,
	std::string cano_PA,
	std::string heading_direction,
	std::string heading_parallel,
	std::string cano_parallel_EE,
	core::Size loop_size,
	core::Size start_res,
	core::Size end_res);

core::Size
WriteToDB_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size new_sheet_id,
	core::Size old_sheet_id);

core::Real
WriteToDB_shortest_dis_between_facing_aro_in_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id,
	core::pose::Pose const & pose,
	utility::vector1<core::Size>	all_distinct_sheet_ids,
	Size	min_num_strands_in_sheet_);

core::Size
WriteToDB_starting_loop(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose & dssp_pose,
	core::Size	sandwich_PK_id_counter,
	core::Size	sw_can_by_sh_id,
	std::string	tag,
	core::Size	max_starting_loop_size_);

core::Size
WriteToDB_sandwich(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	core::Size	sandwich_PK_id_counter,
	std::string	tag,
	core::Size	sw_can_by_sh_id,
	core::Size	sheet_id,
	std::string	sheet_antiparellel,
	core::Size	sandwich_bs_id,
	std::string	strand_is_at_edge,
	core::Size component_size,
	core::Size	residue_begin,
	core::Size	residue_end);


core::Size
WriteToDB_sandwich_by_AA_w_direction(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_w_center_000,
	core::Size	sw_can_by_sh_id,
	core::Size	sheet_id,
	core::Size	residue_begin,
	core::Size	residue_end);

core::Size
WriteToDB_sw_can_by_sh	(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Size sw_can_by_sh_PK_id_counter,
	std::string tag,
	core::Size sw_can_by_sh_id_counter,
	core::Size sheet_id,
	core::Size num_strands_from_sheet);

core::Size
WriteToDB_sw_res_size(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id);


core::Size
WriteToDB_topology_candidate (
	StructureID struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size	sw_can_by_sh_id);

void
WriteToDB_turn_AA(
	core::pose::Pose const & pose,
	core::Size sw_can_by_sh_id,
	core::Size i,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	std::string turn_type);


std::string
WriteToDB_turn_type(
	core::pose::Pose const & pose,
	core::Size sw_can_by_sh_id,
	core::Size start_res,
	core::Size end_res,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::Real	allowed_deviation_for_turn_type_id_);

core::Size
WriteToDB_whether_sw_is_not_connected_with_continuous_atoms(
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session,
	core::Size sw_can_by_sh_id,
	std::string sw_is_not_connected_with_continuous_atoms);


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* WriteToDBFromSandwichFeatures_HH_ */
