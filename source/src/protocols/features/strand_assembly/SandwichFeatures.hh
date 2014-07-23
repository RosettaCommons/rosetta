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
/// @author Doo Nam Kim (doonam.kim@gmail.com, started with Tim Jacobs' code)

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

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
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

	void process_decoy(
		core::pose::Pose & dssp_pose,
		core::scoring::ScoreFunction const&
	) const;


private:

	core::Real
	allowed_deviation_for_turn_type_id_;

	core::Real
	CB_b_factor_cutoff_for_electrostatic_interactions_;

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

	core::Size
	min_num_strands_to_deal_;

	core::Size
	max_num_strands_to_deal_;

	core::Size
	min_res_in_strand_;

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


	core::Size
	max_H_in_extracted_sw_loop_;	//	definition: maximum allowable number of helix residues in extracted sandwich loop

	core::Size
	max_E_in_extracted_sw_loop_;	//	definition: maximum allowable number of E residues in extracted sandwich loop

	core::Real
	max_abs_inter_strand_dihedral_to_not_be_same_direction_strands_;

	core::Real
	max_inter_strand_angle_to_not_be_same_direction_strands_;

	core::Size
	max_starting_loop_size_;

	core::Size
	max_ending_loop_size_;

	core::Size
	max_num_sw_per_pdb_;

	core::Real
	min_N_O_dis_between_two_sheets_;

	core::Real
	min_N_H_O_angle_between_two_sheets_;

	core::Size
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
	write_resfile_to_minimize_too_much_hydrophobic_surface_;

	bool
	write_resfile_to_minimize_too_many_core_heading_FWY_on_core_strands_;

	bool
	write_resfile_to_minimize_too_many_core_heading_FWY_on_edge_strands_;

	bool
	write_strand_AA_distribution_files_;


}; // class SandwichFeatures : public protocols::features::FeaturesReporter

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
