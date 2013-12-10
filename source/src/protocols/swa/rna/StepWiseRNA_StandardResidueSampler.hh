// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/StepWiseRNA_StandardResidueSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_StandardResidueSampler_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_StandardResidueSampler_HH

#include <protocols/swa/rna/StepWiseRNA_StandardResidueSampler.fwd.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <protocols/swa/rna/screener/AtrRepScreener.fwd.hh>
#include <protocols/swa/rna/screener/ChainBreakScreener.fwd.hh>
#include <protocols/swa/rna/screener/ChainClosableScreener.fwd.hh>
#include <protocols/swa/rna/O2PrimePacker.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSelection.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <protocols/rotamer_sampler/RotamerBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_StandardResidueSampler: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseRNA_StandardResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_StandardResidueSampler();

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

	public:

		void
		set_extra_tag( std::string const setting ){ extra_tag_ = setting; }

		void
		set_centroid_screen( bool const setting ){ centroid_screen_ = setting; }

		void
		set_VDW_atr_rep_screen( bool const setting ){ VDW_atr_rep_screen_ = setting; }

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		void set_native_rmsd_screen( bool const & setting ){ native_rmsd_screen_ = setting; }

		void set_native_screen_rmsd_cutoff( core::Real const & setting ){ native_screen_rmsd_cutoff_ = setting; }

		void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

		void set_verbose( bool const & setting ){ verbose_ = setting; }

		void set_perform_o2prime_pack( bool const & setting ){ perform_o2prime_pack_ = setting; }

		void set_cluster_rmsd( core::Real const & setting ){ cluster_rmsd_ = setting; }

		void set_allow_bulge_at_chainbreak( bool const & setting ){ allow_bulge_at_chainbreak_ = setting; }

		void set_pose_list( utility::vector1< PoseOP > &	pose_list );

		utility::vector1< PoseOP > & pose_list();

		core::io::silent::SilentFileDataOP & silent_file_data();

		void output_pose_list( std::string const final_sampler_output_silent_file ) const;

		void set_num_pose_kept( core::Size const & num_pose_kept ){ num_pose_kept_ = num_pose_kept; }

		void
		set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener );

		void
		set_parin_favorite_output( bool const & setting ){ parin_favorite_output_ = setting ; }

		void
		set_allow_syn_pyrimidine(  bool const & setting ){ allow_syn_pyrimidine_ = setting; }

		void
		set_distinguish_pucker( bool const & setting ){ distinguish_pucker_ = setting ; }

		void
		set_use_green_packer( bool const & setting ){ use_green_packer_ = setting ; }

		void
		set_finer_sampling_at_chain_closure( bool const & setting ){ finer_sampling_at_chain_closure_ = setting; }

		void
		set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

		void
		set_reinitialize_CCD_torsions( bool const & setting ){ reinitialize_CCD_torsions_ = setting; }

		void
		set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener );

		void
		set_extra_chi( bool const & setting ){ extra_chi_ = setting; }

		void
		set_extra_epsilon_rotamer( bool const & setting ){ extra_epsilon_rotamer_ = setting; }

		void
		set_extra_beta_rotamer( bool const & setting ){ extra_beta_rotamer_ = setting; }

		void
		set_sample_both_sugar_base_rotamer( bool const & setting ){ sample_both_sugar_base_rotamer_ = setting; }

		void
		set_include_torsion_value_in_tag( bool const & setting ){ include_torsion_value_in_tag_ = setting; }

		void
		set_rebuild_bulge_mode( bool const & setting ){ rebuild_bulge_mode_ = setting; }

		void
		set_combine_long_loop_mode( bool const & setting ){ combine_long_loop_mode_ = setting; }

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
		initialize_poses_and_screeners( core::pose::Pose & pose );

		bool
		break_early_for_integration_tests();

		void
		output_count_data();

		Size
		get_num_pose_kept();

		rotamer_sampler::RotamerBaseOP
		setup_rotamer_sampler( pose::Pose const & pose ) const;

		bool
		bulge_variant_decision( core::pose::Pose & pose, Real const & delta_atr_score );

		void
		apply_bulge_variant( core::pose::Pose & pose ) const;

	private:

		StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010

		Size const moving_res_;
		Size const moving_suite_;
		bool const is_prepend_;
		bool const is_internal_;
		Size const actually_moving_res_;
		Size const gap_size_;
		utility::vector1 < core::Size > const working_moving_partition_pos_;
		Size const num_nucleotides_;
		bool const is_dinucleotide_;
		bool const close_chain_to_distal_;
		Size const five_prime_chain_break_res_;

		Size const last_append_res_;
		Size const last_prepend_res_;
		Real const atom_atom_overlap_dist_cutoff_;
		std::string extra_tag_;

		utility::vector1< PoseOP > pose_list_;
		core::scoring::ScoreFunctionOP scorefxn_;
		StepWiseRNA_CountStruct count_data_;

		std::string silent_file_;
		core::Real const bin_size_; /*ALWAYS 20!!!*/
		core::Size num_pose_kept_;
		core::Real cluster_rmsd_;
		bool verbose_;
		bool native_rmsd_screen_;
		core::Real native_screen_rmsd_cutoff_;

		bool perform_o2prime_pack_;
		bool use_green_packer_;
		bool allow_bulge_at_chainbreak_;
		bool integration_test_mode_;

		screener::AtrRepScreenerOP atr_rep_screener_;
		screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;
		screener::ChainClosableScreenerOP chain_closable_screener_;
		screener::ChainBreakScreenerOP chain_break_screener_;
		screener::StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;
		pose::PoseOP screening_pose_;
		StepWiseRNA_PoseSelectionOP pose_selection_;

		O2PrimePackerOP o2prime_packer_;

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
		bool rebuild_bulge_mode_;
		bool combine_long_loop_mode_;
		bool output_pdb_;
		bool choose_random_;
		Size num_random_samples_;
		bool force_centroid_interaction_;
		bool use_phenix_geo_;
		bool kic_sampling_;

		rotamer_sampler::RotamerBaseOP rotamer_sampler_;

		};

} //rna
} //swa
} //protocols

#endif
