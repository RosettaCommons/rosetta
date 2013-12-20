// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/StepWiseRNA_Modeler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_Modeler_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_Modeler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.fwd.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>
#include <string>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_Modeler: public protocols::moves::Mover  {

	public:

		//constructor
		StepWiseRNA_Modeler( scoring::ScoreFunctionOP scorefxn );

		//constructor
		StepWiseRNA_Modeler( Size const sample_res, scoring::ScoreFunctionOP scorefxn );

		//destructor
		~StepWiseRNA_Modeler();

		StepWiseRNA_Modeler( StepWiseRNA_Modeler const & src );

		StepWiseRNA_ModelerOP clone_modeler() const;

		StepWiseRNA_Modeler & operator=( StepWiseRNA_Modeler const & src );

	public:

		void set_moving_res_and_reset( Size const moving_res );

		virtual void apply( pose::Pose & pose );

		virtual std::string get_name() const;

		void set_job_parameters( StepWiseRNA_JobParametersCOP job_parameters );

		void set_native_pose( pose::PoseCOP );

		pose::PoseCOP get_native_pose();

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_sampler_num_pose_kept( Size const & setting ){ sampler_num_pose_kept_ = setting; }

		void set_num_pose_minimize( Size const & setting ){ num_pose_minimize_ = setting; }

		Size get_num_sampled(){ return num_sampled_; }

		void set_sampler_native_screen_rmsd_cutoff( Real const & setting ){ sampler_native_screen_rmsd_cutoff_ = setting; }

		void set_cluster_rmsd( Real const & setting ){ cluster_rmsd_ = setting; }

		void set_native_edensity_score_cutoff( Real const & setting ){ native_edensity_score_cutoff_ = setting; }

		void set_sampler_native_rmsd_screen( bool const & setting ){ sampler_native_rmsd_screen_ = setting; }

		void set_o2prime_screen( bool const & setting ){ o2prime_screen_ = setting; }

		void set_verbose( bool const & setting ){ verbose_ = setting; }

		void set_distinguish_pucker( bool const & setting ){ distinguish_pucker_ = setting; }

		void set_finer_sampling_at_chain_closure( bool const & setting ){ finer_sampling_at_chain_closure_ = setting; }

		void set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

		void set_allow_syn_pyrimidine( bool const & setting ){ allow_syn_pyrimidine_ = setting; }

		void set_syn_chi_res_list( utility::vector1< core::Size > const & setting ){ syn_chi_res_list_ = setting;}

		void set_extra_chi( bool const & setting ){ extra_chi_ = setting; }

		void set_use_phenix_geo( bool const & setting ){ use_phenix_geo_ = setting; }

		void set_virtual_sugar_legacy_mode( bool const & setting ){ virtual_sugar_legacy_mode_ = setting; }

		void set_kic_sampling( bool const & setting ){ kic_sampling_ = setting; }

		void set_kic_sampling_if_relevant( bool const & setting ){ kic_sampling_if_relevant_ = setting; }

		void set_centroid_screen( bool const & setting ){ centroid_screen_ = setting; }

		void set_VDW_atr_rep_screen( bool const & setting ){ VDW_atr_rep_screen_ = setting; }

		void set_force_centroid_interaction( bool const & setting ){ force_centroid_interaction_ = setting; }

		void set_choose_random( bool const & setting ){ choose_random_ = setting; }

		void set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }

		void set_skip_sampling( bool const & setting ){ skip_sampling_ = setting; }

		void set_perform_minimize( bool const & setting ){ perform_minimize_ = setting; }

		void set_minimize_and_score_sugar( bool const & setting ){ minimize_and_score_sugar_ = setting; }

		void set_minimize_and_score_native_pose( bool const & setting ){ minimize_and_score_native_pose_ = setting; }

		void set_rm_virt_phosphate( bool const & setting ){ rm_virt_phosphate_ = setting; }

		void set_VDW_rep_screen_info( utility::vector1< std::string > const & setting ){ VDW_rep_screen_info_ = setting; }

		void set_VDW_rep_alignment_RMSD_CUTOFF( Real const & setting ){ VDW_rep_alignment_RMSD_CUTOFF_ = setting; }

		void set_output_pdb( bool const & setting ){ output_pdb_ = setting; }

		void set_output_minimized_pose_list( bool const & setting ){ output_minimized_pose_list_ = setting; }

		void set_fixed_res( utility::vector1< Size > const & setting ){ fixed_res_ = setting; minimize_res_.clear(); }

		void set_minimize_res( utility::vector1< Size > const & setting ){ minimize_res_ = setting; fixed_res_.clear(); }

		void set_rmsd_res_list( utility::vector1< Size > const & setting ){ rmsd_res_list_ = setting; }

		// additional options that are not shared with ERRASER (yet)
		void set_VDW_rep_delete_matching_res( utility::vector1< std::string > const & setting ){ VDW_rep_delete_matching_res_ = setting ; }

		void set_VDW_rep_screen_physical_pose_clash_dist_cutoff( bool const & setting ){ VDW_rep_screen_physical_pose_clash_dist_cutoff_ = setting; }

		void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

		void set_allow_bulge_at_chainbreak( bool const & setting ){ allow_bulge_at_chainbreak_ = setting; }

		void set_parin_favorite_output( bool const & setting ){ parin_favorite_output_ = setting; }

		void set_reinitialize_CCD_torsions( bool const & setting ){ reinitialize_CCD_torsions_ = setting; }

		void set_sampler_extra_epsilon_rotamer( bool const & setting ){ sampler_extra_epsilon_rotamer_ = setting; }

		void set_sampler_extra_beta_rotamer( bool const & setting ){ sampler_extra_beta_rotamer_ = setting; }

		void set_sample_both_sugar_base_rotamer( bool const & setting ){ sample_both_sugar_base_rotamer_ = setting; }

		void set_sampler_include_torsion_value_in_tag( bool const & setting ){ sampler_include_torsion_value_in_tag_ = setting; }

		void set_combine_long_loop_mode( bool const & setting ){ combine_long_loop_mode_ = setting; }

		void set_do_not_sample_multiple_virtual_sugar( bool const & setting ){ do_not_sample_multiple_virtual_sugar_ = setting; }

		void set_sample_ONLY_multiple_virtual_sugar( bool const & setting ){ sample_ONLY_multiple_virtual_sugar_ = setting; }

		void set_sampler_assert_no_virt_sugar_sampling( bool const & setting ){ sampler_assert_no_virt_sugar_sampling_ = setting; }

		void set_sampler_try_sugar_instantiation( bool const & setting ){ sampler_try_sugar_instantiation_ = setting; }

		void set_allow_base_pair_only_centroid_screen( bool const & setting ){ allow_base_pair_only_centroid_screen_ = setting; }

		// this is new, not in ERRASER (swa_rna_analytical_closure)
		void set_minimizer_perform_o2prime_pack( bool const & setting ){ minimizer_perform_o2prime_pack_ = setting; }

		void set_minimizer_output_before_o2prime_pack( bool const & setting ){ minimizer_output_before_o2prime_pack_ = setting; }

		void set_minimizer_rename_tag( bool const & setting ){ minimizer_rename_tag_ = setting; }

		void
		set_minimizer_extra_minimize_res( utility::vector1< core::Size > setting ){ minimizer_extra_minimize_res_ = setting; }

		void
		set_minimizer_allow_variable_bond_geometry( bool const setting ){ minimizer_allow_variable_bond_geometry_ = setting; }

		void
		set_minimizer_vary_bond_geometry_frequency( core::Real const setting ) { minimizer_vary_bond_geometry_frequency_ = setting; }


		StepWiseRNA_JobParametersOP
		setup_job_parameters_for_swa_with_full_model_info( utility::vector1< Size > moving_res,
																											 core::pose::Pose & pose );

		StepWiseRNA_JobParametersOP
		setup_job_parameters_for_swa( utility::vector1< Size > moving_res, pose::Pose const & pose );

		void
		output_pose(pose::Pose & pose, std::string const & out_tag, std::string const out_silent_file ) const;

	private:

		void initialize_variables();

		void initialize_job_parameters_and_root( pose::Pose & pose );

		void
		do_residue_sampling( pose::Pose & pose, utility::vector1< PoseOP > & pose_list );

		bool
		sampling_successful( utility::vector1< PoseOP > & pose_list );

		void
		do_minimizing( pose::Pose & pose, utility::vector1< PoseOP > & pose_list );

		void
		add_to_pose_list( utility::vector1< pose::PoseOP > & pose_list, pose::Pose const & pose, std::string const pose_tag ) const;

		void
		setup_root_based_on_full_model_info( pose::Pose & pose, StepWiseRNA_JobParametersCOP & job_parameters );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_;
		pose::PoseCOP native_pose_;

		utility::vector1< Size > moving_res_list_;
		utility::vector1< Size > fixed_res_;
		utility::vector1< Size > minimize_res_;
		utility::vector1< Size > rmsd_res_list_;
		scoring::ScoreFunctionOP scorefxn_;
		std::string silent_file_;
		Size sampler_num_pose_kept_;
		Size num_pose_minimize_;
		Size num_sampled_;
		Real sampler_native_screen_rmsd_cutoff_;
		Real cluster_rmsd_;
		Real native_edensity_score_cutoff_;
		bool sampler_native_rmsd_screen_;
		bool o2prime_screen_;
		bool verbose_;
		bool distinguish_pucker_;
		bool finer_sampling_at_chain_closure_;
		bool PBP_clustering_at_chain_closure_;
		bool allow_syn_pyrimidine_;
		bool extra_chi_;
		bool use_phenix_geo_;
		bool virtual_sugar_legacy_mode_;
		bool kic_sampling_;
		bool kic_sampling_if_relevant_;
		bool centroid_screen_;
		bool VDW_atr_rep_screen_;
		bool force_centroid_interaction_;
		bool choose_random_;
		Size num_random_samples_;
		bool skip_sampling_;
		bool perform_minimize_;
		bool minimize_and_score_sugar_;
		bool minimize_and_score_native_pose_;
		bool rm_virt_phosphate_;

		utility::vector1< std::string > VDW_rep_screen_info_;
		Real VDW_rep_alignment_RMSD_CUTOFF_;
		bool output_pdb_;
		bool output_minimized_pose_list_;

		// additional options that are not shared with ERRASER (yet)
		utility::vector1< std::string > VDW_rep_delete_matching_res_;
		bool VDW_rep_screen_physical_pose_clash_dist_cutoff_;
		bool integration_test_mode_;
		bool allow_bulge_at_chainbreak_;
		bool parin_favorite_output_;
		bool reinitialize_CCD_torsions_;
		bool sampler_extra_epsilon_rotamer_;
		bool sampler_extra_beta_rotamer_;
		bool sample_both_sugar_base_rotamer_;
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
		Real minimizer_vary_bond_geometry_frequency_;
		utility::vector1< core::Size > minimizer_extra_minimize_res_;
		utility::vector1< core::Size > syn_chi_res_list_;

		StepWiseRNA_MinimizerOP stepwise_rna_minimizer_;
		kinematics::MoveMapOP minimize_move_map_;

		screener::StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;
		screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;

		toolbox::AllowInsertOP allow_insert_;
	};

} //rna
} //swa
} //protocols

#endif
