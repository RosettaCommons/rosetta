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

#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSamplerWrapper.fwd.hh>
#include <protocols/swa/rna/SugarModeling.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {

class StepWiseRNA_ResidueSampler: public protocols::moves::Mover {

public:

	//constructor!
	StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

	//destructor -- necessary?
	~StepWiseRNA_ResidueSampler();

	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	void
	set_centroid_screen( bool const setting ){ centroid_screen_ = setting; }

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

	utility::vector1< core::pose::PoseOP > &
	get_pose_list();

	//core::io::silent::SilentFileDataOP & silent_file_data();

	void output_pose_list( std::string const final_sampler_output_silent_file ) const;

	void set_num_pose_kept( core::Size const & num_pose_kept );

	void
	set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener );

	void
	set_parin_favorite_output( bool const & setting ){ parin_favorite_output_ = setting ; }

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
	set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener );

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
	set_combine_long_loop_mode( bool const & setting ){ combine_long_loop_mode_ = setting; }

	void
	set_do_not_sample_multiple_virtual_sugar( bool const & setting ){ do_not_sample_multiple_virtual_sugar_ = setting; }

	void
	set_sample_ONLY_multiple_virtual_sugar( bool const & setting ){ sample_ONLY_multiple_virtual_sugar_ = setting; }

	void
	set_assert_no_virt_sugar_sampling( bool const & setting ){ assert_no_virt_sugar_sampling_ = setting; }

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
	set_virtual_sugar_legacy_mode( bool const & setting ) { virtual_sugar_legacy_mode_ = setting;	}

	void
	set_virtual_sugar_keep_base_fixed( bool const & setting ) { virtual_sugar_keep_base_fixed_ = setting;	}

	void
	set_kic_sampling( bool const & setting ) { kic_sampling_ = setting;	}

	void
	set_try_sugar_instantiation( bool const & setting ) { try_sugar_instantiation_ = setting;	}

private:

	void
	initialize_scorefunctions();

	void
	output_options();

	void
	check_res_not_bulged();

	void
	standard_sampling_WRAPPER( core::pose::Pose & pose );

	void
	floating_base_sampling( core::pose::Pose & pose );

	bool
	instantiate_any_virtual_sugars( core::pose::Pose & pose );

private:

	StepWiseRNA_JobParametersCOP job_parameters_;

	utility::vector1< pose::PoseOP > pose_list_;

	core::scoring::ScoreFunctionOP scorefxn_;

	std::string silent_file_;
	std::string output_filename_;
	core::Size num_pose_kept_;
	core::Real cluster_rmsd_;
	bool verbose_;
	bool native_rmsd_screen_;
	core::Real native_screen_rmsd_cutoff_;

	bool perform_o2prime_pack_;
	bool const use_green_packer_;
	bool allow_bulge_at_chainbreak_;
	bool integration_test_mode_;

	bool parin_favorite_output_;
	bool centroid_screen_;
	bool VDW_atr_rep_screen_;
	bool allow_syn_pyrimidine_;
	bool distinguish_pucker_;
	bool build_pose_from_scratch_;
	bool finer_sampling_at_chain_closure_;
	bool PBP_clustering_at_chain_closure_;
	bool reinitialize_CCD_torsions_;
	bool extra_epsilon_rotamer_;
	bool extra_beta_rotamer_;
	bool extra_chi_;
	bool sample_both_sugar_base_rotamer_;
	bool include_torsion_value_in_tag_;
	bool combine_long_loop_mode_;
	bool do_not_sample_multiple_virtual_sugar_;
	bool sample_ONLY_multiple_virtual_sugar_;
	bool assert_no_virt_sugar_sampling_;
	bool output_pdb_;
	bool choose_random_;
	Size num_random_samples_;
	bool force_centroid_interaction_;
	bool use_phenix_geo_;
	bool virtual_sugar_legacy_mode_;
	bool virtual_sugar_keep_base_fixed_;
	bool kic_sampling_;
	bool try_sugar_instantiation_;

	screener::StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;
	screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;

	StepWiseRNA_VirtualSugarSamplerWrapperOP virtual_sugar_sampler_wrapper_;

};

}
} //swa
} // protocols

#endif
