// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/WriteToFileFromSandwichFeatures.hh
/// @brief Write to a file after SandwichFeatures
/// @author Doo Nam Kim (doonam.kim@gmail.com)
/// @overview

#ifndef INCLUDED_protocols_features_strand_assembly_WriteToFileFromSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_WriteToFileFromSandwichFeatures_HH

//Devel
#include <protocols/features/strand_assembly/CheckForSandwichFeatures.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>
#include <protocols/features/strand_assembly/WriteToDBFromSandwichFeatures.hh>
//#include <protocols/features/strand_assembly/SandwichFeatures.hh>

using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {

	utility::vector1<Size>
	get_vector_of_loop_AA_distribution(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		string loop_kind);

	core::Size
	prepare_and_write_number_of_electrostatic_interactions_of_residues_to_files	(
		std::string tag,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
		core::Real	distance_cutoff_for_electrostatic_interactions_,
		core::Real	CB_b_factor_cutoff_for_electrostatic_interactions_,
		core::Size	min_primary_seq_distance_diff_for_electrostatic_interactions_,
		bool	write_electrostatic_interactions_of_surface_residues_in_a_strand_,
		bool	write_electrostatic_interactions_of_all_residues_in_a_strand_,
		bool	write_electrostatic_interactions_of_all_residues_,
		core::Size	rkde_PK_id_counter,
		core::Size	rkde_in_strands_PK_id_counter);

	core::Size
	write_AA_distribution_with_direction_to_a_file(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session);


	core::Size
	write_AA_distribution_without_direction_to_a_file(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session);

	core::Size
	write_AA_kind_to_a_file(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size	sw_can_by_sh_id,
		core::Size	sw_res_size);

	core::Size
	write_beta_sheet_capping_info_to_a_file(
		string	tag,
		core::pose::Pose const & pose,
		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
		int	primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping_,
		int	primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping_,
		int	primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping_,
		int	primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping_);

	core::Size
	write_chain_B_resNum_to_a_file(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::Size sw_can_by_sh_id	// was 'vec_sw_can_by_sh_id[ii]'
		);

	core::Size
	write_heading_direction_of_all_AA_in_a_strand_to_a_file(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh);

	core::Size
	write_number_of_electrostatic_interactions_of_residues_to_files	(
		std::string tag,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose,
		std::string	dssp_code,
		std::string	heading_direction,
		core::Real	distance_cutoff_for_electrostatic_interactions_,
		core::Real	CB_b_factor_cutoff_for_electrostatic_interactions_,
		core::Size	min_primary_seq_distance_diff_for_electrostatic_interactions_);

	core::Size
	write_p_aa_pp_of_AAs_to_a_file(
		std::string tag,
		core::pose::Pose & dssp_pose);


	core::Size
	write_phi_psi_of_each_residue_to_a_file(
		string	tag,
		core::pose::Pose & dssp_pose,
		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
		bool	write_phi_psi_of_E_,
		bool	write_phi_psi_of_all_,
		core::Size	max_num_sw_per_pdb_,
		StructureID struct_id,
		utility::sql_database::sessionOP  db_session,
		core::Real min_CA_CA_dis_,
		core::Real max_CA_CA_dis_);


	core::Size
	write_rama_of_AAs_to_a_file(
		std::string tag,
		core::pose::Pose & dssp_pose);


	core::Size
	write_resfile_to_a_file(
		string  tag,
		StructureID struct_id,
		utility::sql_database::sessionOP  db_session,
		core::pose::Pose const & pose,
		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
		bool  write_resfile_NOT_FWY_on_surface_);

	core::Size
	write_resfile_to_a_file_when_seq_rec_is_bad(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session,
		core::pose::Pose const & pose,
		utility::vector1<SandwichFragment> bs_of_sw_can_by_sh,
		bool	write_resfile_to_minimize_too_much_hydrophobic_surface_,
		bool	write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_,
		bool	write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_);


} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* WriteToFileFromSandwichFeatures_HH_ */
