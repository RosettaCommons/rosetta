// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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

#include <protocols/features/strand_assembly/CheckForSandwichFeatures.hh>
#include <protocols/features/strand_assembly/SandwichFeatures.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>


using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace core;
using namespace std;

bool
change_sheet_id_if_possible(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Real min_C_O_N_angle_);

Size
delete_this_struct_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);


Size
delete_this_sw_can_by_sh_id_from_sw_by_comp(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

utility::vector1<SandwichFragment>
prepare_WriteToDB_sandwich(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session);

void
WriteToDB_AA_to_terminal_loops (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sandwich_PK_id_counter,
	Size sw_can_by_sh_id,
	string tag,
	bool starting_loop,
	Size residue_begin,
	Size residue_end);


Size
WriteToDB_ending_loop(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sandwich_PK_id_counter,
	Size sw_can_by_sh_id,
	string tag,
	Size max_starting_loop_size_);

Size
WriteToDB_long_strand_id_in_each_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

Size
WriteToDB_avg_b_factor_CB_at_each_component (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id);

Size
WriteToDB_dihedral_angle_between_core_strands_across_facing_sheets (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id);

Size
WriteToDB_dssp_ratio_in_sw(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sw_can_by_sh_id);

Size
WriteToDB_hydrophobic_ratio_net_charge (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);


Size
WriteToDB_min_avg_dis_between_sheets_by_cen_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	core::pose::Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size min_num_strands_in_sheet_);


Size
WriteToDB_min_dis_between_sheets_by_all_res (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	core::pose::Pose & dssp_pose,
	utility::vector1<Size> all_distinct_sheet_ids);

Size
WriteToDB_number_of_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	Size max_num_sw_per_pdb_,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_);

void
WriteToDB_number_of_core_heading_W_in_sw(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

Size
WriteToDB_number_of_core_heading_LWY_in_core_strands_in_sw(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

void
WriteToDB_number_of_core_heading_FWY_in_sw(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);


Size
WriteToDB_number_of_core_heading_aro_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size current_bs_id,
	Size closest_bs_id);

Size
WriteToDB_number_of_core_heading_charged_AAs_in_a_pair_of_edge_strands (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Size current_bs_id,
	Size closest_bs_id);

Size
WriteToDB_number_of_edge_strands_in_each_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

Size
WriteToDB_number_of_sheets_that_surround_this_sheet(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id,
	Size num_of_sheets_that_surround_this_sheet);

Size
WriteToDB_number_strands_in_each_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

Size
WriteToDB_prolines_that_seem_to_prevent_aggregation(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	Real wt_for_pro_in_starting_loop_,
	Real wt_for_pro_in_1st_inter_sheet_loop_,
	Real wt_for_pro_in_3rd_inter_sheet_loop_);

Size
WriteToDB_ratio_of_core_heading_FWY_in_sw(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	core::pose::Pose const & pose);

Size
WriteToDB_rkde(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size rkde_PK_id_counter,
	string tag,
	Size residue_number,
	string residue_type);

Size
WriteToDB_rkde_in_strands(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size rkde_in_strands_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size residue_number,
	string residue_type,
	string heading_direction);

Size
WriteToDB_sheet (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_PK_id_counter,
	Size sheet_id,
	Size segment_id);

Size
WriteToDB_sheet_antiparallel(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sheet_id,
	string antiparallel);


Size
WriteToDB_sheet_connectivity(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sandwich_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	string loop_kind,
	Size intra_sheet_con_id,
	Size inter_sheet_con_id,
	string LR,
	string cano_LR,
	string PA_by_preceding_E,
	string PA_by_following_E,
	string cano_PA,
	string heading_direction,
	string heading_parallel,
	string cano_parallel_EE,
	Size loop_size,
	Size start_res,
	Size end_res);

Size
WriteToDB_sheet_id(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size new_sheet_id,
	Size old_sheet_id);

Real
WriteToDB_shortest_dis_between_facing_aro_in_sw (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	core::pose::Pose const & pose,
	utility::vector1<Size> all_distinct_sheet_ids,
	Size min_num_strands_in_sheet_);

Size
WriteToDB_starting_loop(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose & dssp_pose,
	Size sandwich_PK_id_counter,
	Size sw_can_by_sh_id,
	string tag,
	Size max_starting_loop_size_);

Size
Run_WriteToDB_sandwich(
	string tag,
	core::pose::Pose & dssp_pose,
	utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
	Size max_num_sw_per_pdb_,
	StructureID struct_id,
	utility::sql_database::sessionOP  db_session,
	Real min_CA_CA_dis_,
	Real max_CA_CA_dis_,
	Size sandwich_PK_id_counter);


Size
WriteToDB_sandwich(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	Size sandwich_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id,
	Size sheet_id,
	string sheet_antiparellel,
	Size sandwich_bs_id,
	string strand_is_at_edge,
	Size component_size,
	Size residue_begin,
	Size residue_end);


Size
WriteToDB_sandwich_by_AA_w_direction(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	core::pose::Pose const & pose,
	core::pose::Pose const & pose_w_center_000,
	Size sw_can_by_sh_id,
	Size sheet_id,
	Size residue_begin,
	Size residue_end);

Size
WriteToDB_sw_can_by_sh (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_PK_id_counter,
	string tag,
	Size sw_can_by_sh_id_counter,
	Size sheet_id,
	Size num_strands_from_sheet);

Size
WriteToDB_sw_res_size(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);


Size
WriteToDB_topology_candidate (
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id);

void
WriteToDB_turn_AA(
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id,
	Size i,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	string turn_type);


string
WriteToDB_turn_type(
	core::pose::Pose const & pose,
	Size sw_can_by_sh_id,
	Size start_res,
	Size end_res,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Real allowed_deviation_for_turn_type_id_);

Size
WriteToDB_whether_sw_is_not_connected_with_continuous_atoms(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session,
	Size sw_can_by_sh_id,
	string sw_is_not_connected_with_continuous_atoms);


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* WriteToDBFromSandwichFeatures_HH_ */
