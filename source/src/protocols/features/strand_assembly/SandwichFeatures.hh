// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file SandwichFeatures.hh
/// @brief
/// @author Doo Nam Kim (started with Tim Jacobs' code)

#ifndef INCLUDED_protocols_features_strand_assembly_SandwichFeatures_hh
#define INCLUDED_protocols_features_strand_assembly_SandwichFeatures_hh

//Unit
#include <protocols/features/strand_assembly/SandwichFeatures.fwd.hh>

//External
#include <boost/uuid/uuid.hpp>

//Protocols
#include <protocols/features/FeaturesReporter.hh>

//Devel
#include <protocols/features/strand_assembly/SandwichFragment.hh>

//Utility
#include <utility/vector1.hh>

// for string return
#include <string>

// for get_sw_can_by_sh_id, get_central_residues_in_each_of_two_edge_strands
#include <vector>


// for parse_my_tag
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

class SandwichFeatures : public protocols::features::FeaturesReporter
{

public:

	SandwichFeatures();

	virtual
	std::string
	type_name() const
	{
		return "SandwichFeatures";
	}

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	void
	parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;
	
	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
		core::pose::Pose const & pose, //core::pose::Pose & pose, // dropped 'const' for dssp info addition
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	utility::vector1<SandwichFragment>
	get_full_strands(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	bool
	check_whether_strand_i_is_in_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size segment_id);

	utility::vector1<SandwichFragment>
	get_current_strands_in_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	utility::vector1<SandwichFragment>
	get_all_strands_in_sheet_i(
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

	utility::vector1<core::Size>
	get_chain_B_resNum(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size	sw_can_by_sh_id);

	core::Size
	get_num_of_sheets_that_surround_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size	sheet_id);


	std::string
	get_tag(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	core::Size
	get_num_strands_in_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);

	utility::vector1<SandwichFragment>
	get_full_strands_from_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);

	utility::vector1<core::Size>
	get_distinct_sheet_id_from_sheet_table(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	core::Size
	find_sheet	(
		core::pose::Pose const & pose,
		SandwichFragment strand_i,
		SandwichFragment strand_j,
		bool antiparalell // if false, try to find a sheet in parallel way
		);

	std::vector<core::Real>
	cal_dis_angle_to_find_sheet( // calculate distance and angle to find sheet
		core::pose::Pose const & pose,
		core::Size res_i_0,
		core::Size res_i_1,
		core::Size res_i_2,
		core::Size res_j_0,
		core::Size res_j_1,
		core::Size res_j_2);

	bool
	see_whether_sheets_can_be_combined(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size i_sheet,
		core::Size j_sheet);

	core::Size
	update_sheet_id(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size new_sheet_id,
		core::Size old_sheet_id);

	void
	update_num_of_sheets_that_surround_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id,
		core::Size num_of_sheets_that_surround_this_sheet);

	void
	update_sheet_antiparallel(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id,
		std::string antiparallel);

	std::string
	get_sheet_antiparallel_info(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);

	bool
	change_sheet_id_if_possible(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose);

	core::Real
	check_sw_by_dis(
		core::pose::Pose const & pose,
		SandwichFragment strand_i,
		SandwichFragment strand_j,
		bool antiparalell);

	std::string
	see_whether_sheet_is_antiparallel(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size i_sheet);

	bool
	check_strand_too_closeness	(
		core::pose::Pose const & pose,
		SandwichFragment strand_i,
		SandwichFragment strand_j);

	core::Real
	get_avg_dis_strands(
		core::pose::Pose const & pose,
		SandwichFragment strand_i,
		SandwichFragment strand_j);

	core::Real
	get_closest_distance_between_strands(
		core::pose::Pose const & pose,
		SandwichFragment strand_i,
		SandwichFragment strand_j);

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

	core::Size round(
		core::Real x);

	core::Real round_to_real(
		core::Real x);

	core::Real calculate_dihedral_w_4_resnums(
		core::pose::Pose const & pose,
		core::Size res1_sheet_i,		core::Size res2_sheet_i,		core::Size res1_sheet_j,		core::Size res2_sheet_j);



	// See whether this strand is an edge strand without 'sheet_antiparallel' info
	std::string
	is_this_strand_at_edge	(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id,
		core::Size residue_begin,
		core::Size residue_end);

	std::string
	is_this_strand_at_edge_by_looking_db(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size residue_begin);


	bool
	check_whether_this_sheet_is_too_short(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_i);

	std::pair<core::Size, core::Size>
	get_central_residues_in_each_of_two_edge_strands(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size sheet_i);

	core::Real
	get_shortest_among_4_vals(
		core::Real arr_dis_inter_sheet[]);

	int
	judge_facing(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size sheet_i,
		core::Size sheet_j);

	core::Size
	write_to_sheet (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_PK_id_counter,
		core::Size sheet_id,
		core::Size segment_id);

	core::Size
	write_to_sw_can_by_sh	(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_PK_id_counter,
		std::string tag,
		core::Size sw_can_by_sh_id_counter,
		core::Size sheet_id,
		core::Size num_strands_from_sheet);

	utility::vector1<SandwichFragment>
	prepare_to_fill_sw_by_components(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);



	std::vector<Size>
	get_cen_res_in_other_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sw_can_by_sh_id,
		core::Size sheet_id);

	std::vector<Size>
	get_cen_residues_in_this_sheet(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sheet_id);

	std::vector<core::Size>
	count_AA(
		core::pose::Pose const & pose,
		core::Size residue_begin,
		core::Size residue_end);

	std::vector<core::Size>
	count_AA_w_direction(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		core::Size sw_can_by_sh_id,
		core::Size sheet_id,
		core::Size residue_begin,
		core::Size residue_end);

	core::Size
	fill_sw_by_components(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		core::Size	sw_by_components_PK_id_counter,
		std::string	tag,
		core::Size	sw_can_by_sh_id,
		core::Size	sheet_id,
		std::string	sheet_antiparellel,
		core::Size	sw_by_components_bs_id,
		std::string	strand_is_at_edge,
		core::Size component_size,
		core::Size	residue_begin,
		core::Size	residue_end);

	core::Size
	update_sw_by_components_by_AA_w_direction(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		core::Size	sw_can_by_sh_id,
		core::Size	sheet_id,
		core::Size	residue_begin,
		core::Size	residue_end);
	
	utility::vector1<Size>
	get_vec_sw_can_by_sh_id(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	utility::vector1<Size>
	get_vec_distinct_sheet_id(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size	sw_can_by_sh_id);


	core::Size
	get_size_sw_by_components_PK_id(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id);

	std::pair<core::Size, core::Size>
	get_starting_res_for_connecting_strands(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id,
		core::Size former_res_end);

	std::pair<core::Size, core::Size>
	get_next_starting_res_for_connecting_strands(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id,
		core::Size former_ending_res);

	core::Size
	update_sheet_connectivity(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size sw_by_components_PK_id_counter,
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
	delete_this_sw_can_by_sh_id_from_sw_by_comp(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id);
		
	core::Size
	get_segment_id(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size all_strands_index);
	
	core::Size
	get_num_of_distinct_sheet_id(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	core::Real
	absolute_vec (numeric::xyzVector<core::Real> vector);

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

	std::string
	check_heading_direction ( // positive,,
		core::pose::Pose & dssp_pose,
		core::Size preceding_E,
		core::Size following_E,
		std::string check_heading_direction_by_);

	bool
	check_helix_existence(
	core::pose::Pose const & pose);

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

	bool
	check_whether_same_direction_strands_connect_two_sheets_or_a_loop(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size start_res,
		core::Size next_start_res);

	bool
	check_whether_hairpin_connects_short_strand(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size start_res,
		core::Size next_start_res);

	void
	add_AA_to_terminal_loops (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size	sw_by_components_PK_id_counter,
		core::Size	sw_can_by_sh_id,
		std::string tag,
		bool starting_loop,
		core::Size residue_begin,
		core::Size residue_end);

	core::Size
	add_starting_loop(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size	sw_by_components_PK_id_counter,
		core::Size	sw_can_by_sh_id,
		std::string	tag);

	core::Size
	add_ending_loop(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size	sw_by_components_PK_id_counter,
		core::Size	sw_can_by_sh_id,
		std::string	tag);

	core::Size
	add_dssp_ratio_in_sw(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size sw_can_by_sh_id);

	void
	add_number_of_inward_pointing_W_in_sw(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sw_can_by_sh_id);

	void
	add_number_of_inward_pointing_LWY_in_core_strands_in_sw(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sw_can_by_sh_id);
		
	core::Size
	add_sw_res_size(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sw_can_by_sh_id);

	core::Size
	mark_sw_which_is_not_connected_with_continuous_atoms(
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sw_can_by_sh_id,
		std::string sw_is_not_connected_with_continuous_atoms);


	core::Size
	add_num_strands_in_each_sw (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size	sw_can_by_sh_id);

	core::Size
	add_num_edge_strands_in_each_sw (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size	sw_can_by_sh_id);

	core::Size
	report_hydrophobic_ratio_net_charge (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size	sw_can_by_sh_id);

	core::Size
	report_dihedral_angle_between_core_strands_across_facing_sheets (
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		core::Size	sw_can_by_sh_id);


	bool
	check_whether_this_pdb_should_be_excluded (
		std::string tag); // I don't know how to correctly extract beta-sandwich from 1W8N for now

	core::Size
	cal_num_of_sheets_that_surround_this_sheet (
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		utility::vector1<core::Size>	all_distinct_sheet_ids,
		core::Size sheet_id);

	core::Real
	cal_min_dis_between_sheets_by_cen_res (
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose & dssp_pose,
		core::Size sheet_id_1,
		core::Size sheet_id_2);

	std::pair<core::Size, core::Size>	
	get_current_bs_id_and_closest_edge_bs_id_in_different_sheet (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		core::Size sw_can_by_sh_id,
		core::Size sheet_id,
		core::Size residue_begin,
		core::Size residue_end);

	core::Size
	report_number_of_inward_pointing_charged_AAs_in_a_pair_of_edge_strands (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id,
		core::Size current_bs_id,
		core::Size closest_bs_id);

	core::Size
	report_number_of_inward_pointing_aro_AAs_in_a_pair_of_edge_strands (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id,
		core::Size current_bs_id,
		core::Size closest_bs_id);

	core::Size
	identify_sheet_id_by_residue_end(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size residue_end);

	// used for judge_facing
	utility::vector1<SandwichFragment>
	get_start_end_res_num_in_the_longest_strand(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sheet_id);


	bool
	check_whether_sheets_are_connected_with_near_bb_atoms(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose & dssp_pose,
		core::Size sw_can_by_sh_id);

	std::string
	check_whether_sw_is_not_connected_with_continuous_atoms(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose & dssp_pose,
		core::Size sw_can_by_sh_id);


	utility::vector1<core::Size>
	get_vector_AA_distribution_w_direction (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string heading_direction, // like core_heading, surface_heading
		std::string strand_location // like edge_strand, core_strand
		);


	utility::vector1<core::Size>
	get_vector_AA_distribution_wo_direction (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string loop_kind
		);


	utility::vector1<core::Size>
	get_vec_AA_kind (
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id);

		
	bool
	check_whether_sw_by_sh_id_still_alive(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::Size sw_can_by_sh_id);


	void
	report_turn_type(
		core::pose::Pose const & pose,
		core::Size sw_can_by_sh_id,
		core::Size start_res,
		core::Size end_res,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);
	
 
private:

	core::Size
	min_num_strands_to_deal_;
	
	core::Size
	max_num_strands_to_deal_;
		
	core::Size
	min_res_in_strand_;

	core::Real
	min_CA_CA_dis_;

	core::Real
	max_CA_CA_dis_;

	core::Real
	min_O_N_dis_;

	core::Real
	max_O_N_dis_;

	core::Real
	min_C_O_N_angle_;

	core::Real
	min_sheet_dis_;

	core::Real
	max_sheet_dis_;

	core::Real
	max_sheet_angle_with_cen_res_in_smaller_sheet_and_two_terminal_res_in_larger_sheet_;

	core::Real
	min_sheet_angle_by_four_term_cen_res_;

	core::Real
	max_sheet_angle_by_four_term_cen_res_;

	core::Real
	min_sheet_torsion_cen_res_;

	core::Real
	max_sheet_torsion_cen_res_;

	core::Size
	min_num_strands_in_sheet_; //  definition: a sheet with < 3 strands will be ignored

	core::Real
	min_inter_sheet_dis_CA_CA_;

	core::Real
	max_inter_sheet_dis_CA_CA_;

	bool
	extract_sandwich_;
		
	bool	
	write_chain_B_resnum_;

	bool
	no_helix_in_pdb_;

	core::Size
	max_H_in_extracted_sw_loop_;	//	definition: maximum allowable number of helix residues in extracted sandwich loop

	core::Size
	max_E_in_extracted_sw_loop_;	//	definition: maximum allowable number of E residues in extracted sandwich loop

	bool
	exclude_sandwich_that_is_linked_w_same_direction_strand_;

	core::Real
	max_inter_strand_angle_to_not_be_same_direction_strands_;

	core::Real
	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_;

	bool	
	write_phi_psi_of_all_;

	bool
	write_phi_psi_of_E_;

	core::Size
	max_starting_loop_size_;

	core::Size
	max_ending_loop_size_;


	core::Size
	max_num_sw_per_pdb_;

	std::string
	check_N_to_C_direction_by_; // 1) PE: preceding_E's CA-CB vector,
								// 2) FE: following_E's CA-CB vector,
								// 3) CBs: preceding_E's CB to following_E's CB vector

	bool
	do_not_connect_sheets_by_loops_;
	
	core::Real	
	check_canonicalness_cutoff_;

	bool
	count_AA_with_direction_;

	core::Real
	inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_;

	bool
	exclude_desinated_pdbs_;

	bool
	exclude_sandwich_that_has_near_backbone_atoms_between_sheets_;

	core::Real
	min_N_O_dis_between_two_sheets_;

	core::Real
	min_N_H_O_angle_between_two_sheets_;

	bool
	write_all_info_files_;

	bool
	write_AA_kind_files_;

	bool
	write_AA_distribution_files_w_direction_;

	bool
	write_AA_distribution_files_wo_direction_;

}; // class SandwichFeatures : public protocols::features::FeaturesReporter

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
