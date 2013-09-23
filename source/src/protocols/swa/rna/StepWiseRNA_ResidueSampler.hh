// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software res and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Parin Sripakdeevong
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_ResidueSampler_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_ResidueSampler_HH

//#include <numeric/xyzMatrix.hh>
//#include <numeric/xyzVector.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualRiboseSampler.hh>
#include <protocols/swa/rna/FloatingBaseChainClosureJobParameter.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBaseSamplerUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.hh>     //Feb 02, 2012: Need this to pass rosetta_tools/python_cc_reader/test_all_headers_compile.py
#include <protocols/swa/rna/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>
#include <protocols/rotamer_sampler/RotamerBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>


namespace protocols {
namespace swa {
namespace rna {

class StepWiseRNA_ResidueSampler: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

	//destructor -- necessary?
	~StepWiseRNA_ResidueSampler();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	void
	set_centroid_screen( bool const setting ){ centroid_screen_ = setting; }

	void
	set_allow_base_pair_only_centroid_screen( bool const setting ){ allow_base_pair_only_centroid_screen_ = setting; }

	void
	set_VDW_atr_rep_screen( bool const setting ){ VDW_atr_rep_screen_ = setting; }

	void
	set_silent_file( std::string const & setting );

	void
	set_output_filename( std::string const & output_filename );

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	void
	set_native_rmsd_screen( bool const & setting );

	void
	set_native_screen_rmsd_cutoff( core::Real const & setting );

	void
	set_integration_test_mode( bool const & setting );

	void
	set_verbose( bool const & setting );

	void
	set_perform_o2prime_pack( bool const & setting );

	void
	set_cluster_rmsd( core::Real const & setting );

	void
	set_allow_bulge_at_chainbreak( bool const & setting );

	utility::vector1< PoseOP > &
	get_pose_data_list();

	core::io::silent::SilentFileDataOP & silent_file_data();

	void output_pose_data_list( std::string const final_sampler_output_silent_file ) const;

	void set_num_pose_kept( core::Size const & num_pose_kept );

	void
	set_base_centroid_screener( StepWiseRNA_BaseCentroidScreenerOP & screener );

	void
	set_parin_favorite_output( bool const & setting ){ parin_favorite_output_ = setting ; }

	void
	set_floating_base( bool const & setting ){ floating_base_ = setting ; }

	void
	set_allow_syn_pyrimidine(  bool const & setting ){ allow_syn_pyrimidine_ = setting; }

	void
	set_distinguish_pucker( bool const & setting ){ distinguish_pucker_ = setting ; }

	void
	set_finer_sampling_at_chain_closure( bool const & setting ){ finer_sampling_at_chain_closure_ = setting; }

	void
	set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

	void
	set_reinitialize_CCD_torsions( bool const & setting ){ reinitialize_CCD_torsions_ = setting; }

	void
	set_user_input_VDW_bin_screener( StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener ){ user_input_VDW_bin_screener_ = user_input_VDW_bin_screener; }

	void
	set_extra_epsilon_rotamer( bool const & setting ){ extra_epsilon_rotamer_ = setting; }

	void
	set_extra_beta_rotamer( bool const & setting ){ extra_beta_rotamer_ = setting; }

	void
	set_extra_chi( bool const & setting ){ extra_chi_ = setting; }

	void
	set_sample_both_sugar_base_rotamer( bool const & setting ){ sample_both_sugar_base_rotamer_ = setting; }

	void
	set_include_torsion_value_in_tag( bool const & setting ){ include_torsion_value_in_tag_ = setting; }

	void
	set_rebuild_bulge_mode( bool const & setting ){ rebuild_bulge_mode_ = setting; }

	void
	set_combine_long_loop_mode( bool const & setting ){ combine_long_loop_mode_ = setting; }

	void
	set_do_not_sample_multiple_virtual_sugar( bool const & setting ){ do_not_sample_multiple_virtual_sugar_ = setting; }

	void
	set_sample_ONLY_multiple_virtual_sugar( bool const & setting ){ sample_ONLY_multiple_virtual_sugar_ = setting; }

	void
	set_assert_no_virt_ribose_sampling( bool const & setting ){ assert_no_virt_ribose_sampling_ = setting; }

	void
	set_output_pdb( bool const setting ){ output_pdb_ = setting; }

	void
	set_choose_random( bool const & setting ) {	choose_random_ = setting;	}

	void
	set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }

	void
	set_force_centroid_interaction ( bool const & setting ) {		force_centroid_interaction_ = setting;	}

	void
	set_use_phenix_geo( bool const & setting ) { use_phenix_geo_ = setting;	}

	void
	set_kic_sampling( bool const & setting ) { kic_sampling_ = setting;	}


private:

	void
	initialize_scorefunctions();

	void
	Copy_CCD_torsions( core::pose::Pose & pose, core::pose::Pose const & template_pose ) const;

	void
	Copy_CCD_torsions_general( core::pose::Pose & pose, core::pose::Pose const & template_pose, core::Size const five_prime_res, core::Size const three_prime_res ) const;

	bool
	check_loop_closed( core::pose::Pose const & pose );

	bool
	Chain_break_screening( core::pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & constraint_scorefxn );

	bool
	Chain_break_screening_general( core::pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn, core::Size const five_prime_res );

	void
	standard_sampling_WRAPPER( core::pose::Pose & pose,
														FloatingBaseChainClosureJobParameter const & prev_sugar_FB_JP,
														FloatingBaseChainClosureJobParameter const & curr_sugar_FB_JP,
														FloatingBaseChainClosureJobParameter const & five_prime_CB_sugar_FB_JP,
														FloatingBaseChainClosureJobParameter const & three_prime_CB_sugar_FB_JP );

	void
	standard_sampling( core::pose::Pose & pose, utility::vector1< PoseOP > & pose_data_list, std::string const sugar_tag );


	void
	floating_base_sampling( core::pose::Pose & pose, FloatingBaseChainClosureJobParameter const & prev_sugar_FB_JP );

	void
	get_base_atr_rep_score( core::pose::Pose const & pose, core::Real & base_atr_score, core::Real & base_rep_score );

	bool
	Full_atom_van_der_Waals_screening_REPLICATE(
																		core::pose::Pose & current_pose_screen,
																		core::Real const & base_rep_score,
																		core::Real const &  base_atr_score,
																		core::Real & delta_rep_score,
																		core::Real & delta_atr_score,
																		core::Size const & gap_size,
																		bool const & Is_internal );


	bool
	Full_atom_van_der_Waals_screening(
																		core::pose::Pose & current_pose_screen,
																		core::Real const & base_rep_score,
																		core::Real const &  base_atr_score,
																		core::Real & delta_rep_score,
																		core::Real & delta_atr_score,
																		core::Size const & gap_size,
																		bool const & Is_internal );

	void
	initialize_o2prime_packer_task( core::pose::Pose const & pose );

	void
	initialize_o2prime_green_packer( core::pose::Pose & pose );

	void
	sample_o2prime_hydrogen( core::pose::Pose & pose, core::pose::Pose & pose_with_original_HO2prime_torsion );

	core::Real
	Pose_selection_by_full_score( utility::vector1< PoseOP > & pose_data_list, core::pose::Pose & current_pose, std::string const & tag );

	bool
	apply_bulge_variant( core::pose::Pose & pose, core::Real const & delta_atr_score );

	void
	Update_pose_data_list( std::string const & tag, utility::vector1< PoseOP > & pose_data_list, core::pose::Pose const & current_pose, core::Real const & current_score ) const;

	void
	cluster_pose_data_list( utility::vector1< PoseOP > & pose_data_list );

	std::string
	create_tag( std::string const & prestring, core::Size const i ) const;

	core::kinematics::Stub
	get_reference_stub( core::Size const reference_res, core::pose::Pose const & pose ) const;


	std::string //silly function to convert to real to string
	create_torsion_value_string( core::Real const & torsion_value ) const;

	std::string
	create_rotamer_string( core::pose::Pose const & pose ) const;

	/////////////////////////////////////function related to sampling/setup virtual sugar /////////////////////////////////////////////////////////////////
	bool
	Is_previous_sugar_virtual( core::pose::Pose const & pose ) const;

	bool
	Is_current_sugar_virtual( core::pose::Pose const & pose ) const;

	bool
	Is_five_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) const;

	bool
	Is_three_prime_chain_break_sugar_virtual( core::pose::Pose const & pose ) const;

	utility::vector1< PoseOP >
	previous_floating_base_chain_closure(
		core::pose::Pose & viewer_pose,
		FloatingBaseChainClosureJobParameter const & FB_job_params,
		std::string const name
	);

	void output_count_data();

	rotamer_sampler::RotamerBaseOP
	setup_rotamer_sampler( core::pose::Pose const & pose ) const;

private:

	StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010

	core::io::silent::SilentFileDataOP sfd_;
	utility::vector1< PoseOP > pose_data_list_;

	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;
	core::scoring::ScoreFunctionOP chainbreak_scorefxn_;
	core::scoring::ScoreFunctionOP sampling_scorefxn_;
	core::scoring::ScoreFunctionOP o2prime_pack_scorefxn_;

	utility::vector1 < core::Size > working_rmsd_res_;

	SillyCountStruct count_data_;
	std::string silent_file_;
	std::string output_filename_;
	core::Real const bin_size_; /*ALWAYS 20!!!*/
	core::Real const rep_cutoff_;
	core::Size num_pose_kept_;
	core::Size const multiplier_;
	core::Real cluster_rmsd_;
	bool verbose_;
	bool native_rmsd_screen_;
	core::Real native_screen_rmsd_cutoff_;

	bool perform_o2prime_pack_;
	core::pack::task::PackerTaskOP o2prime_pack_task_;
	protocols::simple_moves::GreenPackerOP o2prime_green_packer_;
	bool const use_green_packer_;
	bool allow_bulge_at_chainbreak_;
	bool integration_test_mode_;
	bool floating_base_;

	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;

	bool parin_favorite_output_;
	bool centroid_screen_;
	bool allow_base_pair_only_centroid_screen_;
	bool VDW_atr_rep_screen_;
	bool allow_syn_pyrimidine_;
	bool distinguish_pucker_;
	bool build_pose_from_scratch_;
	core::Real current_score_cutoff_;
	bool finer_sampling_at_chain_closure_;
	bool PBP_clustering_at_chain_closure_;
	bool reinitialize_CCD_torsions_;
	bool extra_epsilon_rotamer_;
	bool extra_beta_rotamer_;
	bool extra_chi_;
	bool sample_both_sugar_base_rotamer_;
	bool include_torsion_value_in_tag_;
	bool rebuild_bulge_mode_;
	bool combine_long_loop_mode_;
	bool do_not_sample_multiple_virtual_sugar_;
	bool sample_ONLY_multiple_virtual_sugar_;
	bool assert_no_virt_ribose_sampling_;
	bool output_pdb_;
	bool choose_random_;
	Size num_random_samples_;
	bool force_centroid_interaction_;
	bool use_phenix_geo_;
	bool kic_sampling_;

	StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;
};

}
} //swa
} // protocols

#endif
