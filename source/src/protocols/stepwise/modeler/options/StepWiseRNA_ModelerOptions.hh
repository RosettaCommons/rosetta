// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/options/StepWiseRNA_ModelerOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_StepWiseRNA_ModelerOptions_HH
#define INCLUDED_protocols_stepwise_modeler_rna_StepWiseRNA_ModelerOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

	// multiple inheritance -- bad form -- but will replace with composition later, perhaps.
	class StepWiseRNA_ModelerOptions: public virtual basic::resource_manager::ResourceOptions {

	public:

		//constructor
		StepWiseRNA_ModelerOptions();

		StepWiseRNA_ModelerOptions( StepWiseRNA_ModelerOptions const & src );

		//destructor
		~StepWiseRNA_ModelerOptions();

	public:

		StepWiseRNA_ModelerOptionsOP clone() const;

		/// @brief Describe this instance to a given output stream
		virtual
		void
		show( std::ostream & ) const{}

		/// @brief Initialize from the recursive "tag" structure.
		virtual
		void
		parse_my_tag( utility::tag::TagCOP ){}

		/// @brief The class name (its type) for a particular ResourceOptions instance.
		/// This function allows for better error message delivery.
		virtual
		std::string
		type() const{ return "StepWiseRNA_ModelerOptions";}

	public:

		void
		initialize_from_command_line();

		core::Real const & native_edensity_score_cutoff() const { return native_edensity_score_cutoff_; }
		void set_native_edensity_score_cutoff( core::Real const & setting ){ native_edensity_score_cutoff_ = setting; }

		bool const & o2prime_legacy_mode() const { return o2prime_legacy_mode_; }
		void set_o2prime_legacy_mode( bool const & setting ){ o2prime_legacy_mode_ = setting; }

		bool const & allow_virtual_o2prime_hydrogens() const { return allow_virtual_o2prime_hydrogens_; }
		void set_allow_virtual_o2prime_hydrogens( bool const & setting ){ allow_virtual_o2prime_hydrogens_ = setting; }

		bool const & sampler_perform_phosphate_pack() const { return sampler_perform_phosphate_pack_; }
		void set_sampler_perform_phosphate_pack( bool const & setting ){ sampler_perform_phosphate_pack_ = setting; }

		bool const & distinguish_pucker() const { return distinguish_pucker_; }
		void set_distinguish_pucker( bool const & setting ){ distinguish_pucker_ = setting; }

		bool const & finer_sampling_at_chain_closure() const { return finer_sampling_at_chain_closure_; }
		void set_finer_sampling_at_chain_closure( bool const & setting ){ finer_sampling_at_chain_closure_ = setting; }

		bool const & PBP_clustering_at_chain_closure() const { return PBP_clustering_at_chain_closure_; }
		void set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

		bool const & allow_syn_pyrimidine() const { return allow_syn_pyrimidine_; }
		void set_allow_syn_pyrimidine( bool const & setting ){ allow_syn_pyrimidine_ = setting; }

		bool const & extra_chi() const { return extra_chi_; }
		void set_extra_chi( bool const & setting ){ extra_chi_ = setting; }

		bool const & use_phenix_geo() const { return use_phenix_geo_; }
		void set_use_phenix_geo( bool const & setting ){ use_phenix_geo_ = setting; }

		bool const & virtual_sugar_legacy_mode() const { return virtual_sugar_legacy_mode_; }
		void set_virtual_sugar_legacy_mode( bool const & setting ){ virtual_sugar_legacy_mode_ = setting; }

		bool const & virtual_sugar_keep_base_fixed() const { return virtual_sugar_keep_base_fixed_; }
		void set_virtual_sugar_keep_base_fixed( bool const & setting ){ virtual_sugar_keep_base_fixed_ = setting; }

		bool const & virtual_sugar_do_minimize() const { return virtual_sugar_do_minimize_; }
		void set_virtual_sugar_do_minimize( bool const & setting ){ virtual_sugar_do_minimize_ = setting; }

		bool const & kic_modeler_if_relevant() const { return kic_modeler_if_relevant_; }
		void set_kic_modeler_if_relevant( bool const & setting ){ kic_modeler_if_relevant_ = setting; }

		bool const & force_centroid_interaction() const { return force_centroid_interaction_; }
		void set_force_centroid_interaction( bool const & setting ){ force_centroid_interaction_ = setting; }

		bool const & minimize_and_score_sugar() const { return minimize_and_score_sugar_; }
		void set_minimize_and_score_sugar( bool const & setting ){ minimize_and_score_sugar_ = setting; }

		bool const & minimize_and_score_native_pose() const { return minimize_and_score_native_pose_; }
		void set_minimize_and_score_native_pose( bool const & setting ){ minimize_and_score_native_pose_ = setting; }

		bool const & rm_virt_phosphate() const { return rm_virt_phosphate_; }
		void set_rm_virt_phosphate( bool const & setting ){ rm_virt_phosphate_ = setting; }

		core::Real const & VDW_rep_alignment_RMSD_CUTOFF() const { return VDW_rep_alignment_RMSD_CUTOFF_; }
		void set_VDW_rep_alignment_RMSD_CUTOFF( core::Real const & setting ){ VDW_rep_alignment_RMSD_CUTOFF_ = setting; }

		core::Real const & sampler_max_centroid_distance() const { return sampler_max_centroid_distance_; }
		void set_sampler_max_centroid_distance( core::Real const & setting ){ sampler_max_centroid_distance_ = setting; }

		utility::vector1< std::string > const & VDW_rep_delete_matching_res() const { return VDW_rep_delete_matching_res_; }
		void set_VDW_rep_delete_matching_res( utility::vector1< std::string > const & setting ){ VDW_rep_delete_matching_res_ = setting; }

		bool const & VDW_rep_screen_physical_pose_clash_dist_cutoff() const { return VDW_rep_screen_physical_pose_clash_dist_cutoff_; }
		void set_VDW_rep_screen_physical_pose_clash_dist_cutoff( bool const & setting ){ VDW_rep_screen_physical_pose_clash_dist_cutoff_ = setting; }

		bool const & integration_test_mode() const { return integration_test_mode_; }
		void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

		bool const & allow_bulge_at_chainbreak() const { return allow_bulge_at_chainbreak_; }
		void set_allow_bulge_at_chainbreak( bool const & setting ){ allow_bulge_at_chainbreak_ = setting; }

		bool const & parin_favorite_output() const { return parin_favorite_output_; }
		void set_parin_favorite_output( bool const & setting ){ parin_favorite_output_ = setting; }

		bool const & reinitialize_CCD_torsions() const { return reinitialize_CCD_torsions_; }
		void set_reinitialize_CCD_torsions( bool const & setting ){ reinitialize_CCD_torsions_ = setting; }

		bool const & sampler_extra_epsilon_rotamer() const { return sampler_extra_epsilon_rotamer_; }
		void set_sampler_extra_epsilon_rotamer( bool const & setting ){ sampler_extra_epsilon_rotamer_ = setting; }

		bool const & sampler_extra_beta_rotamer() const { return sampler_extra_beta_rotamer_; }
		void set_sampler_extra_beta_rotamer( bool const & setting ){ sampler_extra_beta_rotamer_ = setting; }

		bool const & sampler_include_torsion_value_in_tag() const { return sampler_include_torsion_value_in_tag_; }
		void set_sampler_include_torsion_value_in_tag( bool const & setting ){ sampler_include_torsion_value_in_tag_ = setting; }

		bool const & combine_long_loop_mode() const { return combine_long_loop_mode_; }
		void set_combine_long_loop_mode( bool const & setting ){ combine_long_loop_mode_ = setting; }

		bool const & do_not_sample_multiple_virtual_sugar() const { return do_not_sample_multiple_virtual_sugar_; }
		void set_do_not_sample_multiple_virtual_sugar( bool const & setting ){ do_not_sample_multiple_virtual_sugar_ = setting; }

		bool const & sample_ONLY_multiple_virtual_sugar() const { return sample_ONLY_multiple_virtual_sugar_; }
		void set_sample_ONLY_multiple_virtual_sugar( bool const & setting ){ sample_ONLY_multiple_virtual_sugar_ = setting; }

		bool const & sampler_assert_no_virt_sugar_sampling() const { return sampler_assert_no_virt_sugar_sampling_; }
		void set_sampler_assert_no_virt_sugar_sampling( bool const & setting ){ sampler_assert_no_virt_sugar_sampling_ = setting; }

		bool const & sampler_try_sugar_instantiation() const { return sampler_try_sugar_instantiation_; }
		void set_sampler_try_sugar_instantiation( bool const & setting ){ sampler_try_sugar_instantiation_ = setting; }

		bool const & allow_base_pair_only_centroid_screen() const { return allow_base_pair_only_centroid_screen_; }
		void set_allow_base_pair_only_centroid_screen( bool const & setting ){ allow_base_pair_only_centroid_screen_ = setting; }

		bool const & minimizer_rename_tag() const { return minimizer_rename_tag_; }
		void set_minimizer_rename_tag( bool const & setting ){ minimizer_rename_tag_ = setting; }

		bool const & tether_jump() const { return tether_jump_; }
		void set_tether_jump( bool const & setting ){ tether_jump_ = setting; }

		bool const & turn_off_rna_chem_map_during_optimize() const { return turn_off_rna_chem_map_during_optimize_; }
		void set_turn_off_rna_chem_map_during_optimize( bool const & setting ){ turn_off_rna_chem_map_during_optimize_ = setting; }

	protected:

		void
		initialize_variables();

	private:

		core::Real native_edensity_score_cutoff_;
		bool o2prime_legacy_mode_;
		bool allow_virtual_o2prime_hydrogens_;
		bool sampler_perform_phosphate_pack_;
		bool distinguish_pucker_;
		bool finer_sampling_at_chain_closure_;
		bool PBP_clustering_at_chain_closure_;
		bool allow_syn_pyrimidine_;
		bool extra_chi_;
		bool use_phenix_geo_;
		bool virtual_sugar_legacy_mode_;
		bool virtual_sugar_keep_base_fixed_;
		bool virtual_sugar_do_minimize_;
		bool kic_modeler_if_relevant_;
		bool force_centroid_interaction_;
		bool minimize_and_score_sugar_;
		bool minimize_and_score_native_pose_;
		bool rm_virt_phosphate_;
		core::Real VDW_rep_alignment_RMSD_CUTOFF_;
		core::Real sampler_max_centroid_distance_;
		utility::vector1< std::string > VDW_rep_delete_matching_res_;
		bool VDW_rep_screen_physical_pose_clash_dist_cutoff_;
		bool integration_test_mode_;
		bool allow_bulge_at_chainbreak_;
		bool parin_favorite_output_;
		bool reinitialize_CCD_torsions_;
		bool sampler_extra_epsilon_rotamer_;
		bool sampler_extra_beta_rotamer_;
		bool sampler_include_torsion_value_in_tag_;
		bool combine_long_loop_mode_;
		bool do_not_sample_multiple_virtual_sugar_;
		bool sample_ONLY_multiple_virtual_sugar_;
		bool sampler_assert_no_virt_sugar_sampling_;
		bool sampler_try_sugar_instantiation_;
		bool allow_base_pair_only_centroid_screen_;
		bool minimizer_rename_tag_;
		bool tether_jump_;
		bool turn_off_rna_chem_map_during_optimize_;
	};

} //options
} //modeler
} //stepwise
} //protocols

#endif
