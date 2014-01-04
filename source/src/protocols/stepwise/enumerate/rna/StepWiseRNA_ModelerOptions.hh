// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_enumerate_rna_StepWiseRNA_ModelerOptions_HH
#define INCLUDED_protocols_stepwise_enumerate_rna_StepWiseRNA_ModelerOptions_HH

#include <basic/resource_manager/ResourceOptions.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {

	class StepWiseRNA_ModelerOptions: public basic::resource_manager::ResourceOptions {

	public:

		//constructor
		StepWiseRNA_ModelerOptions();

		StepWiseRNA_ModelerOptions( StepWiseRNA_ModelerOptions const & src );

		//destructor
		~StepWiseRNA_ModelerOptions();

		StepWiseRNA_ModelerOptionsOP clone() const;

	public:

		StepWiseRNA_ModelerOptions &
		operator = ( StepWiseRNA_ModelerOptions const & src );

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
		type() const{ return "StepWiseRNA_ModelOptions";}

	public:

		void
		initialize_from_command_line();

		void
		setup_options_for_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener ) const;

		std::string const & silent_file() const { return silent_file_; }
		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		core::Size const & sampler_num_pose_kept() const { return sampler_num_pose_kept_; }
		void set_sampler_num_pose_kept( core::Size const & setting ){ sampler_num_pose_kept_ = setting; }

		core::Size const & num_pose_minimize() const { return num_pose_minimize_; }
		void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

		bool const & sampler_native_rmsd_screen() const { return sampler_native_rmsd_screen_; }
		void set_sampler_native_rmsd_screen( bool const & setting ){ sampler_native_rmsd_screen_ = setting; }

		core::Real const & sampler_native_screen_rmsd_cutoff() const { return sampler_native_screen_rmsd_cutoff_; }
		void set_sampler_native_screen_rmsd_cutoff( core::Real const & setting ){ sampler_native_screen_rmsd_cutoff_ = setting; }

		core::Real const & cluster_rmsd() const { return cluster_rmsd_; }
		void set_cluster_rmsd( core::Real const & setting ){ cluster_rmsd_ = setting; }

		core::Real const & native_edensity_score_cutoff() const { return native_edensity_score_cutoff_; }
		void set_native_edensity_score_cutoff( core::Real const & setting ){ native_edensity_score_cutoff_ = setting; }

		bool const & sampler_perform_o2prime_pack() const { return sampler_perform_o2prime_pack_; }
		void set_sampler_perform_o2prime_pack( bool const & setting ){ sampler_perform_o2prime_pack_ = setting; }

		bool const & use_green_packer() const { return use_green_packer_; }
		void set_use_green_packer( bool const & setting ){ use_green_packer_ = setting; }

		bool const & verbose() const { return verbose_; }
		void set_verbose( bool const & setting ){ verbose_ = setting; }

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

		bool const & kic_sampling_if_relevant() const { return kic_sampling_if_relevant_; }
		void set_kic_sampling_if_relevant( bool const & setting ){ kic_sampling_if_relevant_ = setting; }

		bool const & VDW_atr_rep_screen() const { return VDW_atr_rep_screen_; }
		void set_VDW_atr_rep_screen( bool const & setting ){ VDW_atr_rep_screen_ = setting; }

		bool const & force_centroid_interaction() const { return force_centroid_interaction_; }
		void set_force_centroid_interaction( bool const & setting ){ force_centroid_interaction_ = setting; }

		bool const & choose_random() const { return choose_random_; }
		void set_choose_random( bool const & setting ){ choose_random_ = setting; }

		Size const & num_random_samples() const { return num_random_samples_; };
		void set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }

		bool const & skip_sampling() const { return skip_sampling_; }
		void set_skip_sampling( bool const & setting ){ skip_sampling_ = setting; }

		bool const & perform_minimize() const { return perform_minimize_; }
		void set_perform_minimize( bool const & setting ){ perform_minimize_ = setting; }

		bool const & minimize_and_score_sugar() const { return minimize_and_score_sugar_; }
		void set_minimize_and_score_sugar( bool const & setting ){ minimize_and_score_sugar_ = setting; }

		bool const & minimize_and_score_native_pose() const { return minimize_and_score_native_pose_; }
		void set_minimize_and_score_native_pose( bool const & setting ){ minimize_and_score_native_pose_ = setting; }

		bool const & rm_virt_phosphate() const { return rm_virt_phosphate_; }
		void set_rm_virt_phosphate( bool const & setting ){ rm_virt_phosphate_ = setting; }

		core::Real const & VDW_rep_alignment_RMSD_CUTOFF() const { return VDW_rep_alignment_RMSD_CUTOFF_; }
		void set_VDW_rep_alignment_RMSD_CUTOFF( core::Real const & setting ){ VDW_rep_alignment_RMSD_CUTOFF_ = setting; }

		bool const & output_pdb() const { return output_pdb_; }
		void set_output_pdb( bool const & setting ){ output_pdb_ = setting; }

		bool const & output_minimized_pose_list() const { return output_minimized_pose_list_; }
		void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

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

		bool const & minimizer_perform_o2prime_pack() const { return minimizer_perform_o2prime_pack_; }
		void set_minimizer_perform_o2prime_pack( bool const & setting ){ minimizer_perform_o2prime_pack_ = setting; }

		bool const & minimizer_output_before_o2prime_pack() const { return minimizer_output_before_o2prime_pack_; }
		void set_minimizer_output_before_o2prime_pack( bool const & setting ){ minimizer_output_before_o2prime_pack_ = setting; }

		bool const & minimizer_rename_tag() const { return minimizer_rename_tag_; }
		void set_minimizer_rename_tag( bool const & setting ){ minimizer_rename_tag_ = setting; }

		bool const & minimizer_allow_variable_bond_geometry() const { return minimizer_allow_variable_bond_geometry_; }
		void set_minimizer_allow_variable_bond_geometry( bool const & setting ){ minimizer_allow_variable_bond_geometry_ = setting; }

		core::Real const & minimizer_vary_bond_geometry_frequency() const { return minimizer_vary_bond_geometry_frequency_; }
		void set_minimizer_vary_bond_geometry_frequency( core::Real const & setting ){ minimizer_vary_bond_geometry_frequency_ = setting; }


	private:

		void
		initialize_variables();

		std::string silent_file_;
		core::Size sampler_num_pose_kept_;
		core::Size num_pose_minimize_;
		bool sampler_native_rmsd_screen_;
		core::Real sampler_native_screen_rmsd_cutoff_;
		core::Real cluster_rmsd_;
		core::Real native_edensity_score_cutoff_;
		bool sampler_perform_o2prime_pack_;
		bool use_green_packer_;
		bool verbose_;
		bool distinguish_pucker_;
		bool finer_sampling_at_chain_closure_;
		bool PBP_clustering_at_chain_closure_;
		bool allow_syn_pyrimidine_;
		bool extra_chi_;
		bool use_phenix_geo_;
		bool virtual_sugar_legacy_mode_;
		bool virtual_sugar_keep_base_fixed_;
		bool kic_sampling_if_relevant_;
		bool VDW_atr_rep_screen_;
		bool force_centroid_interaction_;
		bool choose_random_;
		Size num_random_samples_;
		bool skip_sampling_;
		bool perform_minimize_;
		bool minimize_and_score_sugar_;
		bool minimize_and_score_native_pose_;
		bool rm_virt_phosphate_;
		bool output_pdb_;
		bool output_minimized_pose_list_;
		core::Real VDW_rep_alignment_RMSD_CUTOFF_;
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
		bool minimizer_perform_o2prime_pack_;
		bool minimizer_output_before_o2prime_pack_;
		bool minimizer_rename_tag_;
		bool minimizer_allow_variable_bond_geometry_;
		core::Real minimizer_vary_bond_geometry_frequency_;

	};

} //rna
} //enumerate
} //stepwise
} //protocols

#endif
