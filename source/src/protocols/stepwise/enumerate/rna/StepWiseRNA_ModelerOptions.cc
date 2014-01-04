// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/Tracer.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "protocols.stepwise.enumerate.rna.StepWiseRNA_ModelerOptions" );

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {

	/////////////////////////////////////////////////////////////////////////////////////
	//Constructor
	StepWiseRNA_ModelerOptions::StepWiseRNA_ModelerOptions()
	{
		initialize_variables();
	}

	/////////////////////////////////////////////////////////////////////////////////////
	//Destructor
	StepWiseRNA_ModelerOptions::~StepWiseRNA_ModelerOptions()
	{
	}

	/// @brief copy constructor
	StepWiseRNA_ModelerOptions::StepWiseRNA_ModelerOptions( StepWiseRNA_ModelerOptions const & src ) :
		ResourceOptions ( src )
	{
		*this = src;
	}


	/////////////////////////////////////////////////////////////////////////////////////
	// If you add a variable, initialize it here, and include in operator = definition below!
	void
	StepWiseRNA_ModelerOptions::initialize_variables(){
		silent_file_ = "";
		sampler_num_pose_kept_ = 108;
		num_pose_minimize_ = 99999;
		sampler_native_rmsd_screen_ = false;
		sampler_native_screen_rmsd_cutoff_ = 2.0;
		cluster_rmsd_ = 0.5;
		native_edensity_score_cutoff_ = -1;
		sampler_perform_o2prime_pack_ = true;
		use_green_packer_ = false;
		verbose_ = false;
		distinguish_pucker_ = true;
		finer_sampling_at_chain_closure_ = false;
		PBP_clustering_at_chain_closure_ = false;
		allow_syn_pyrimidine_ = false;
		extra_chi_ = false;
		use_phenix_geo_ = false;
		virtual_sugar_legacy_mode_ = false;
		virtual_sugar_keep_base_fixed_ = true;
		kic_sampling_if_relevant_ = false;
		VDW_atr_rep_screen_ = true;
		force_centroid_interaction_ = false;
		choose_random_ = false;
		num_random_samples_ = 1;
		skip_sampling_ = false;
		perform_minimize_ = true;
		minimize_and_score_sugar_ = true;
		minimize_and_score_native_pose_ = false;
		rm_virt_phosphate_ = false;
		VDW_rep_alignment_RMSD_CUTOFF_ = 0.001;
		output_pdb_ = false;
		output_minimized_pose_list_ = false;
		VDW_rep_screen_physical_pose_clash_dist_cutoff_ = false;
		integration_test_mode_ = false;
		allow_bulge_at_chainbreak_ = false;
		parin_favorite_output_ = false;
		reinitialize_CCD_torsions_ = false;
		sampler_extra_epsilon_rotamer_ = false;
		sampler_extra_beta_rotamer_ = false;
		sampler_include_torsion_value_in_tag_ = false;
		combine_long_loop_mode_ = false;
		do_not_sample_multiple_virtual_sugar_ = false;
		sample_ONLY_multiple_virtual_sugar_ = false;
		sampler_assert_no_virt_sugar_sampling_ = false;
		sampler_try_sugar_instantiation_ = false;
		allow_base_pair_only_centroid_screen_ = false;
		minimizer_perform_o2prime_pack_ = true;
		minimizer_output_before_o2prime_pack_ = false;
		minimizer_rename_tag_ = false;
		minimizer_allow_variable_bond_geometry_ = false;
		minimizer_vary_bond_geometry_frequency_ = 0.0;
	}

	/// @brief clone the options
	StepWiseRNA_ModelerOptionsOP
	StepWiseRNA_ModelerOptions::clone() const
	{
		return new StepWiseRNA_ModelerOptions( *this );
	}

	/////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_ModelerOptions &
	StepWiseRNA_ModelerOptions::operator = ( StepWiseRNA_ModelerOptions const & src )
	{
		silent_file_ = src.silent_file_;
		sampler_num_pose_kept_ = src.sampler_num_pose_kept_;
		num_pose_minimize_ = src.num_pose_minimize_;
		sampler_native_rmsd_screen_ = src.sampler_native_rmsd_screen_;
		sampler_native_screen_rmsd_cutoff_ = src.sampler_native_screen_rmsd_cutoff_;
		cluster_rmsd_ = src.cluster_rmsd_;
		native_edensity_score_cutoff_ = src.native_edensity_score_cutoff_;
		sampler_perform_o2prime_pack_ = src.sampler_perform_o2prime_pack_;
		use_green_packer_ = src.use_green_packer_;
		verbose_ = src.verbose_;
		distinguish_pucker_ = src.distinguish_pucker_;
		finer_sampling_at_chain_closure_ = src.finer_sampling_at_chain_closure_;
		PBP_clustering_at_chain_closure_ = src.PBP_clustering_at_chain_closure_;
		allow_syn_pyrimidine_ = src.allow_syn_pyrimidine_;
		extra_chi_ = src.extra_chi_;
		use_phenix_geo_ = src.use_phenix_geo_;
		virtual_sugar_legacy_mode_ = src.virtual_sugar_legacy_mode_;
		virtual_sugar_keep_base_fixed_ = src.virtual_sugar_keep_base_fixed_;
		kic_sampling_if_relevant_ = src.kic_sampling_if_relevant_;
		VDW_atr_rep_screen_ = src.VDW_atr_rep_screen_;
		force_centroid_interaction_ = src.force_centroid_interaction_;
		choose_random_ = src.choose_random_;
		num_random_samples_ = src.num_random_samples_;
		skip_sampling_ = src.skip_sampling_;
		perform_minimize_ = src.perform_minimize_;
		minimize_and_score_sugar_ = src.minimize_and_score_sugar_;
		minimize_and_score_native_pose_ = src.minimize_and_score_native_pose_;
		rm_virt_phosphate_ = src.rm_virt_phosphate_;
		VDW_rep_alignment_RMSD_CUTOFF_ = src.VDW_rep_alignment_RMSD_CUTOFF_;
		output_pdb_ = src.output_pdb_;
		output_minimized_pose_list_ = src.output_minimized_pose_list_;
		VDW_rep_screen_physical_pose_clash_dist_cutoff_ = src.VDW_rep_screen_physical_pose_clash_dist_cutoff_;
		integration_test_mode_ = src.integration_test_mode_;
		allow_bulge_at_chainbreak_ = src.allow_bulge_at_chainbreak_;
		parin_favorite_output_ = src.parin_favorite_output_;
		reinitialize_CCD_torsions_ = src.reinitialize_CCD_torsions_;
		sampler_extra_epsilon_rotamer_ = src.sampler_extra_epsilon_rotamer_;
		sampler_extra_beta_rotamer_ = src.sampler_extra_beta_rotamer_;
		sampler_include_torsion_value_in_tag_ = src.sampler_include_torsion_value_in_tag_;
		combine_long_loop_mode_ = src.combine_long_loop_mode_;
		do_not_sample_multiple_virtual_sugar_ = src.do_not_sample_multiple_virtual_sugar_;
		sample_ONLY_multiple_virtual_sugar_ = src.sample_ONLY_multiple_virtual_sugar_;
		sampler_assert_no_virt_sugar_sampling_ = src.sampler_assert_no_virt_sugar_sampling_;
		sampler_try_sugar_instantiation_ = src.sampler_try_sugar_instantiation_;
		allow_base_pair_only_centroid_screen_ = src.allow_base_pair_only_centroid_screen_;
		minimizer_perform_o2prime_pack_ = src.minimizer_perform_o2prime_pack_;
		minimizer_output_before_o2prime_pack_ = src.minimizer_output_before_o2prime_pack_;
		minimizer_rename_tag_ = src.minimizer_rename_tag_;
		minimizer_vary_bond_geometry_frequency_ = src.minimizer_vary_bond_geometry_frequency_;
		minimizer_allow_variable_bond_geometry_ = src.minimizer_allow_variable_bond_geometry_;

		return *this;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ModelerOptions::initialize_from_command_line(){

		using namespace basic::options::OptionKeys::stepwise::rna;

		sampler_num_pose_kept_ = option[ basic::options::OptionKeys::stepwise::rna::sampler_num_pose_kept ]();
		if ( option[ basic::options::OptionKeys::stepwise::rna::choose_random ]() ) num_pose_minimize_ = 1;
		if ( option[ basic::options::OptionKeys::stepwise::rna::num_pose_minimize ].user() ) num_pose_minimize_ = option[ basic::options::OptionKeys::stepwise::rna::num_pose_minimize ]();

		sampler_native_rmsd_screen_ = option[ OptionKeys::stepwise::rna::sampler_native_rmsd_screen ]();
		sampler_native_screen_rmsd_cutoff_ = option[ OptionKeys::stepwise::rna::sampler_native_screen_rmsd_cutoff ]();
		cluster_rmsd_ =	option[ OptionKeys::stepwise::rna::sampler_cluster_rmsd ]() ;
		native_edensity_score_cutoff_ = option[ OptionKeys::stepwise::rna::native_edensity_score_cutoff]();
		sampler_perform_o2prime_pack_ = option[ OptionKeys::stepwise::rna::sampler_perform_o2prime_pack ]();
		use_green_packer_ = option[ OptionKeys::stepwise::rna::sampler_use_green_packer ]();
		verbose_ = option[ VERBOSE ]();
		distinguish_pucker_ = option[ OptionKeys::stepwise::rna::distinguish_pucker]();
		finer_sampling_at_chain_closure_ = option[ OptionKeys::stepwise::rna::finer_sampling_at_chain_closure]();
		PBP_clustering_at_chain_closure_ = option[ OptionKeys::stepwise::rna::PBP_clustering_at_chain_closure]();
		allow_syn_pyrimidine_ = option[ OptionKeys::stepwise::rna::sampler_allow_syn_pyrimidine ]();
		extra_chi_ = option[ OptionKeys::stepwise::rna::sampler_extra_chi_rotamer]();
		use_phenix_geo_ = option[ basic::options::OptionKeys::rna::corrected_geo ]();
		virtual_sugar_legacy_mode_ = option[ OptionKeys::stepwise::rna::virtual_sugar_legacy_mode ];
		virtual_sugar_keep_base_fixed_ = option[ OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed ]();
		kic_sampling_if_relevant_ = option[ OptionKeys::stepwise::rna::erraser ]();
		VDW_atr_rep_screen_ = option[ OptionKeys::stepwise::rna::VDW_atr_rep_screen ]();
		force_centroid_interaction_ = option[ OptionKeys::stepwise::rna::force_centroid_interaction ]();
		choose_random_ = option[ basic::options::OptionKeys::stepwise::rna::choose_random ]() ;
		num_random_samples_ = option[ basic::options::OptionKeys::stepwise::rna::num_random_samples ]();
		skip_sampling_ = option[ OptionKeys::stepwise::rna::skip_sampling ]();
		perform_minimize_ = option[ OptionKeys::stepwise::rna::minimizer_perform_minimize ]();
		minimize_and_score_sugar_ = option[ basic::options::OptionKeys::stepwise::rna::minimize_and_score_sugar ]();
		minimize_and_score_native_pose_ = option[ OptionKeys::stepwise::rna::minimize_and_score_native_pose ]();
		rm_virt_phosphate_ = option[ OptionKeys::stepwise::rna::rm_virt_phosphate]();
		VDW_rep_alignment_RMSD_CUTOFF_ = option[ OptionKeys::stepwise::rna::VDW_rep_alignment_RMSD_CUTOFF ]();
				VDW_rep_screen_physical_pose_clash_dist_cutoff_ = option[ OptionKeys::stepwise::rna::VDW_rep_screen_physical_pose_clash_dist_cutoff ]();
		// newer options
		VDW_rep_delete_matching_res_ = option[ OptionKeys::stepwise::rna::VDW_rep_delete_matching_res ](); // wait where is this?
		integration_test_mode_ = option[ OptionKeys::stepwise::rna::integration_test ](); //Should set after setting sampler_native_screen_rmsd_cutoff, fast, medium_fast options.
		allow_bulge_at_chainbreak_ = option[ OptionKeys::stepwise::rna::allow_bulge_at_chainbreak ]();
		parin_favorite_output_ = option[ OptionKeys::stepwise::rna::parin_favorite_output ]();
		reinitialize_CCD_torsions_ = option[ OptionKeys::stepwise::rna::reinitialize_CCD_torsions]();
		sampler_extra_epsilon_rotamer_ = option[ OptionKeys::stepwise::rna::sampler_extra_epsilon_rotamer]();
		sampler_extra_beta_rotamer_ = option[ OptionKeys::stepwise::rna::sampler_extra_beta_rotamer]();
		sampler_include_torsion_value_in_tag_ = option[ OptionKeys::stepwise::rna::sampler_include_torsion_value_in_tag]();
		combine_long_loop_mode_ = option[ basic::options::OptionKeys::stepwise::rna::combine_long_loop_mode]();
		do_not_sample_multiple_virtual_sugar_ = option[ OptionKeys::stepwise::rna::do_not_sample_multiple_virtual_sugar]();
		sample_ONLY_multiple_virtual_sugar_ = option[ OptionKeys::stepwise::rna::sample_ONLY_multiple_virtual_sugar]();
		sampler_assert_no_virt_sugar_sampling_ = option[ OptionKeys::stepwise::rna::sampler_assert_no_virt_sugar_sampling ]();
		sampler_try_sugar_instantiation_ = option[ OptionKeys::stepwise::rna::sampler_try_sugar_instantiation ]();
		allow_base_pair_only_centroid_screen_ = option[ OptionKeys::stepwise::rna::allow_base_pair_only_centroid_screen ]();
		minimizer_perform_o2prime_pack_ = option[ OptionKeys::stepwise::rna::minimizer_perform_o2prime_pack ]();
		minimizer_output_before_o2prime_pack_ = option[ OptionKeys::stepwise::rna::minimizer_output_before_o2prime_pack ]();
		minimizer_rename_tag_ = option[ OptionKeys::stepwise::rna::minimizer_rename_tag ]();
		//		minimizer_allow_variable_bond_geometry_ = option[ OptionKeys::stepwise::rna::minimizer_allow_variable_bond_geometry ]();
		//		minimizer_vary_bond_geometry_frequency_ = option[ OptionKeys::stepwise::rna::minimizer_vary_bond_geometry_frequency ]();

		if ( integration_test_mode_ ) num_pose_minimize_ = 1;

	}

	////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ModelerOptions::setup_options_for_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener ) const {
		user_input_VDW_bin_screener->set_VDW_rep_alignment_RMSD_CUTOFF ( VDW_rep_alignment_RMSD_CUTOFF_ );
		user_input_VDW_bin_screener->set_VDW_rep_delete_matching_res( VDW_rep_delete_matching_res_ );
		user_input_VDW_bin_screener->set_physical_pose_clash_dist_cutoff( VDW_rep_screen_physical_pose_clash_dist_cutoff_ );
		user_input_VDW_bin_screener->set_output_pdb( output_pdb_ );
	}


} //rna
} //enumerate
} //stepwise
} //protocols
