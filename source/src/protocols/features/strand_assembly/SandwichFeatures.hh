// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/strand_assembly/SandwichFeatures.hh
/// @brief Extract and analyze beta-sandwich features
/// @author Doo Nam Kim (doonam.kim@gmail.com, started from Tim Jacobs' code)
/// @overview
///  @ task 0: Determine whether we deal with given pdb file
///  @ task 1: Identify all beta-strands
///  @ task 2: Identify all beta-sheets with these strands
///  @ task 3: Identify all beta-sandwiches with these sheets
///   @ task 3-1: Merge sheets to each other if possible
///   @ task 3-2: Make beta-sandwiches with sheets that are ideal only
///    @ task 3-2-1: Exclude if this_sheet_is_surrounded_by_more_than_1_other_sheet
///    @ task 3-2-2: Exclude sheets that are too close to each other
///    @ task 3-2-3: Exclude sheets that are too distant to each other
///    @ task 3-2-4: Exclude sheets that do not face each other
///     @ task 3-2-4-1: Exclude sheets that do not face each other by an angle with two terminal residues and one central residue
///     @ task 3-2-4-2: Exclude sheets that do not face each other by an angle with four terminal residues in two edge strands
///   @ task 3-4: Test canonical sandwich test
///    @ task 3-4-1: Canonical sandwiches need to have low number of helix or strand residues in any loop (beta-hairpin-loop or inter-sheet-loop)
///    @ task 3-4-2: Canonical sandwiches need to not have same-direction-strands as connecting two beta-sheets
///    @ task 3-4-3: Canonical sandwiches should not be beta-barrel obviously
///  @ task 4: Write beta-sandwiches that passed canonical tests into database
///    @ task 4-1: Write hairpin_loop and inter-sheet loop
///    @ task 4-2: Write starting_loop and endng_loop
///    @ task 4-3: Write ratio_hydrophobic_philic/net_charge
///    @ task 4-4: Write total size of sandwich
///  @ task 5: Write resfiles automatically


#ifndef INCLUDED_protocols_features_strand_assembly_SandwichFeatures_hh
#define INCLUDED_protocols_features_strand_assembly_SandwichFeatures_hh

//Devel
#include <protocols/features/strand_assembly/SandwichFeatures.fwd.hh>
#include <protocols/features/strand_assembly/SandwichFragment.hh>
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>

namespace protocols {
namespace features {
namespace strand_assembly {

class SandwichFeatures : public protocols::features::FeaturesReporter
{

public:

	SandwichFeatures();
	~SandwichFeatures();

	virtual
	std::string
	type_name() const
	{
		return "SandwichFeatures";
	}

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief collect all the feature data for the pose
	virtual
	Size
	report_features(
		core::pose::Pose const & pose, //core::pose::Pose & pose, // dropped 'const' for dssp info addition
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);


private:

	core::Real
		allowed_deviation_for_turn_type_id_;

	core::Real
		CB_b_factor_cutoff_for_electrostatic_interactions_;

	bool
		chance_of_being_canonical_sw;

	std::string
		check_N_to_C_direction_by_; // 1) PE: preceding_E's CA-CB vector,
	// 2) FE: following_E's CA-CB vector,
	// 3) CBs: preceding_E's CB to following_E's CB vector


	core::Real
		check_canonicalness_cutoff_;

	bool
		count_AA_with_direction_;

	core::Real
		distance_cutoff_for_electrostatic_interactions_;

	bool
		do_not_connect_sheets_by_loops_;

	bool
		do_not_write_resfile_of_sandwich_that_has_non_canonical_LR_;

	bool
		exclude_sandwich_that_is_linked_w_same_direction_strand_;

	bool
		exclude_sandwich_that_has_non_canonical_LR_;

	bool
		exclude_sandwich_that_has_non_canonical_properties_;

	bool
		exclude_sandwich_that_has_non_canonical_shortest_dis_between_facing_aro_in_sw_;

	bool
		exclude_sandwich_with_SS_bond_;

	bool
		exclude_desinated_pdbs_;

	bool
		exclude_sandwich_that_is_suspected_to_have_not_facing_2_sheets_;

	bool
		exclude_sandwich_that_has_near_backbone_atoms_between_sheets_;



	bool
		extract_sandwich_;

	/// @brief create score-functions for centroid and fullatom level
	core::scoring::ScoreFunctionOP
	generate_scorefxn( bool fullatom = false );


	core::Real
		inter_sheet_distance_to_see_whether_a_sheet_is_surrounded_by_other_sheets_;


	core::Real
		min_CA_CA_dis_;

	core::Real
		max_CA_CA_dis_;

	Size
		min_num_strands_to_deal_;

	Size
		max_num_strands_to_deal_;

	Size
		min_res_in_strand_;

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

	Size
		min_num_strands_in_sheet_; //  definition: a sheet with < 3 strands will be ignored

	core::Real
		min_inter_sheet_dis_CA_CA_;

	core::Real
		max_inter_sheet_dis_CA_CA_;


	Size
		max_H_in_extracted_sw_loop_; // definition: maximum allowable number of helix residues in extracted sandwich loop

	Size
		max_E_in_extracted_sw_loop_; // definition: maximum allowable number of E residues in extracted sandwich loop

	core::Real
		max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_;

	core::Real
		max_inter_strand_angle_to_not_be_same_direction_strands_;

	Size
		max_starting_loop_size_;

	Size
		max_ending_loop_size_;

	Size
		max_num_sw_per_pdb_;

	core::Real
		min_N_O_dis_between_two_sheets_;

	core::Real
		min_N_H_O_angle_between_two_sheets_;

	Size
		min_primary_seq_distance_diff_for_electrostatic_interactions_;

	bool
		no_helix_in_pdb_;

	int
		primary_seq_distance_cutoff_for_beta_sheet_capping_before_N_term_capping_;

	int
		primary_seq_distance_cutoff_for_beta_sheet_capping_after_N_term_capping_;

	int
		primary_seq_distance_cutoff_for_beta_sheet_capping_before_C_term_capping_;

	int
		primary_seq_distance_cutoff_for_beta_sheet_capping_after_C_term_capping_;

	bool
		unit_test_pass_identifier;

	core::Real
		wt_for_pro_in_starting_loop_;

	core::Real
		wt_for_pro_in_1st_inter_sheet_loop_;

	core::Real
		wt_for_pro_in_3rd_inter_sheet_loop_;


	bool
		write_all_info_files_;

	bool
		write_AA_kind_files_;

	bool
		write_beta_sheet_capping_info_;

	bool
		write_chain_B_resnum_;

	bool
		write_electrostatic_interactions_of_surface_residues_in_a_strand_;

	bool
		write_electrostatic_interactions_of_all_residues_in_a_strand_;

	bool
		write_electrostatic_interactions_of_all_residues_;

	bool
		write_heading_directions_of_all_AA_in_a_strand_;

	bool
		write_loop_AA_distribution_files_;

	bool
		write_p_aa_pp_files_;

	bool
		write_phi_psi_of_all_;

	bool
		write_phi_psi_of_E_;

	bool
		write_rama_at_AA_to_files_;

	bool
		write_resfile_;

	bool
		write_resfile_NOT_FWY_on_surface_;

	bool
		write_resfile_to_minimize_too_much_hydrophobic_surface_;

	bool
		write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_;

	bool
		write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_;

	bool
		write_resfile_when_seq_rec_is_bad_;

	bool
		write_strand_AA_distribution_files_;


}; // class SandwichFeatures : public protocols::features::FeaturesReporter

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
