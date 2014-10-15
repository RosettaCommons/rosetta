// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/options/StepWiseRNA_ModelerOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/options/StepWiseRNA_ModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <utility/tag/Tag.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.options.StepWiseRNA_ModelerOptions" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

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

		native_edensity_score_cutoff_ = -1;
		o2prime_legacy_mode_ = false;
		allow_virtual_o2prime_hydrogens_ = false;
		sampler_perform_phosphate_pack_ = true;
		distinguish_pucker_ = true;
		finer_sampling_at_chain_closure_ = false;
		PBP_clustering_at_chain_closure_ = false;
		allow_syn_pyrimidine_ = false;
		extra_chi_ = false;
		use_phenix_geo_ = false;
		virtual_sugar_legacy_mode_ = false;
		virtual_sugar_keep_base_fixed_ = true;
		virtual_sugar_do_minimize_ = true;
		kic_modeler_if_relevant_ = false;
		force_centroid_interaction_ = false;
		minimize_and_score_sugar_ = true;
		minimize_and_score_native_pose_ = false;
		rm_virt_phosphate_ = false;
		VDW_rep_alignment_RMSD_CUTOFF_ = 0.001;
		sampler_max_centroid_distance_ = 0.0;
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
		minimizer_rename_tag_ = false;
		minimizer_allow_variable_bond_geometry_ = false;
		minimizer_vary_bond_geometry_frequency_ = 0.0;
		tether_jump_ = true;
		turn_off_rna_chem_map_during_optimize_ = true;
	}

	/// @brief clone the options
	StepWiseRNA_ModelerOptionsOP
	StepWiseRNA_ModelerOptions::clone() const
	{
		return StepWiseRNA_ModelerOptionsOP( new StepWiseRNA_ModelerOptions( *this ) );
	}

	/////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_ModelerOptions &
	StepWiseRNA_ModelerOptions::operator = ( StepWiseRNA_ModelerOptions const & src )
	{
		native_edensity_score_cutoff_ = src.native_edensity_score_cutoff_;
		o2prime_legacy_mode_ = src.o2prime_legacy_mode_;
		allow_virtual_o2prime_hydrogens_ = src.allow_virtual_o2prime_hydrogens_;
		sampler_perform_phosphate_pack_ = src.sampler_perform_phosphate_pack_;
		distinguish_pucker_ = src.distinguish_pucker_;
		finer_sampling_at_chain_closure_ = src.finer_sampling_at_chain_closure_;
		PBP_clustering_at_chain_closure_ = src.PBP_clustering_at_chain_closure_;
		allow_syn_pyrimidine_ = src.allow_syn_pyrimidine_;
		extra_chi_ = src.extra_chi_;
		use_phenix_geo_ = src.use_phenix_geo_;
		virtual_sugar_legacy_mode_ = src.virtual_sugar_legacy_mode_;
		virtual_sugar_keep_base_fixed_ = src.virtual_sugar_keep_base_fixed_;
		virtual_sugar_do_minimize_ = src.virtual_sugar_do_minimize_;
		kic_modeler_if_relevant_ = src.kic_modeler_if_relevant_;
		force_centroid_interaction_ = src.force_centroid_interaction_;
		minimize_and_score_sugar_ = src.minimize_and_score_sugar_;
		minimize_and_score_native_pose_ = src.minimize_and_score_native_pose_;
		rm_virt_phosphate_ = src.rm_virt_phosphate_;
		VDW_rep_screen_info_ = src.VDW_rep_screen_info_;
		VDW_rep_alignment_RMSD_CUTOFF_ = src.VDW_rep_alignment_RMSD_CUTOFF_;
		sampler_max_centroid_distance_ = src.sampler_max_centroid_distance_;
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
		minimizer_vary_bond_geometry_frequency_ = src.minimizer_vary_bond_geometry_frequency_;
		minimizer_allow_variable_bond_geometry_ = src.minimizer_allow_variable_bond_geometry_;
		tether_jump_ = src.tether_jump_;
		turn_off_rna_chem_map_during_optimize_ = src.turn_off_rna_chem_map_during_optimize_;

		return *this;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_ModelerOptions::initialize_from_command_line(){

		using namespace basic::options::OptionKeys::stepwise::rna;

		native_edensity_score_cutoff_ = option[ OptionKeys::stepwise::rna::native_edensity_score_cutoff]();
		o2prime_legacy_mode_ = option[ OptionKeys::stepwise::rna::o2prime_legacy_mode ]();
		allow_virtual_o2prime_hydrogens_ = option[ OptionKeys::stepwise::rna::allow_virtual_o2prime_hydrogens ]();
		sampler_perform_phosphate_pack_ = option[ OptionKeys::stepwise::rna::sampler_perform_phosphate_pack ]();
		distinguish_pucker_ = option[ OptionKeys::stepwise::rna::distinguish_pucker]();
		finer_sampling_at_chain_closure_ = option[ OptionKeys::stepwise::rna::finer_sampling_at_chain_closure]();
		PBP_clustering_at_chain_closure_ = option[ OptionKeys::stepwise::rna::PBP_clustering_at_chain_closure]();
		allow_syn_pyrimidine_ = option[ OptionKeys::stepwise::rna::sampler_allow_syn_pyrimidine ]();
		extra_chi_ = option[ OptionKeys::stepwise::rna::sampler_extra_chi_rotamer]();
		use_phenix_geo_ = option[ basic::options::OptionKeys::rna::corrected_geo ]();
		virtual_sugar_legacy_mode_ = option[ OptionKeys::stepwise::rna::virtual_sugar_legacy_mode ];
		virtual_sugar_keep_base_fixed_ = option[ OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed ]();
		virtual_sugar_do_minimize_ = option[ OptionKeys::stepwise::rna::virtual_sugar_do_minimize ]();
		kic_modeler_if_relevant_ = option[ OptionKeys::stepwise::rna::erraser ]();
		force_centroid_interaction_ = option[ OptionKeys::stepwise::rna::force_centroid_interaction ]();
		minimize_and_score_sugar_ = option[ basic::options::OptionKeys::stepwise::rna::minimize_and_score_sugar ]();
		minimize_and_score_native_pose_ = option[ OptionKeys::stepwise::rna::minimize_and_score_native_pose ]();
		rm_virt_phosphate_ = option[ OptionKeys::stepwise::rna::rm_virt_phosphate]();
		VDW_rep_screen_info_ = option[ OptionKeys::stepwise::rna::VDW_rep_screen_info ]();
		VDW_rep_alignment_RMSD_CUTOFF_ = option[ OptionKeys::stepwise::rna::VDW_rep_alignment_RMSD_CUTOFF ]();
		sampler_max_centroid_distance_ = option[ OptionKeys::stepwise::rna::sampler_max_centroid_distance ]();
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
		minimizer_rename_tag_ = option[ OptionKeys::stepwise::rna::minimizer_rename_tag ]();
		tether_jump_ = option[ OptionKeys::stepwise::rna::tether_jump ]();
		turn_off_rna_chem_map_during_optimize_ = option[ OptionKeys::stepwise::rna::turn_off_rna_chem_map_during_optimize ]();

	}


} //options
} //modeler
} //stepwise
} //protocols
