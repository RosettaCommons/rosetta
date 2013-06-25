// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/ERRASER_Modeler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_ERRASER_Modeler_HH
#define INCLUDED_protocols_swa_rna_ERRASER_Modeler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/swa/rna/ERRASER_Modeler.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <string>

namespace protocols {
namespace swa {
namespace rna {

	class ERRASER_Modeler: public protocols::moves::Mover  {

	public:

	//constructor
		ERRASER_Modeler( core::Size const sample_res, core::scoring::ScoreFunctionOP scorefxn );

		// should also make the following as an overloaded constructor.
		//ERRASER_Modeler( utility::vector1< core::Size > const moving_res,
		//										 core::scoring::ScoreFunctionOP scorefxn );

		//destructor
		~ERRASER_Modeler();

	public:

		virtual void apply( core::pose::Pose & pose );

		virtual std::string get_name() const;

		void set_job_parameters( StepWiseRNA_JobParametersCOP job_parameters );

		void set_native_pose( core::pose::PoseCOP );

		core::pose::PoseCOP get_native_pose();

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_sampler_num_pose_kept( core::Size const & setting ){ sampler_num_pose_kept_ = setting; }

		void set_num_pose_minimize( core::Size const & setting ){ num_pose_minimize_ = setting; }

		void set_nstruct( core::Size const & setting ){ nstruct_ = setting; }

		core::Size get_num_sampled(){ return num_sampled_; }

		void set_sampler_native_screen_rmsd_cutoff( core::Real const & setting ){ sampler_native_screen_rmsd_cutoff_ = setting; }

		void set_cluster_rmsd( core::Real const & setting ){ cluster_rmsd_ = setting; }

		void set_native_edensity_score_cutoff( core::Real const & setting ){ native_edensity_score_cutoff_ = setting; }

		void set_fast( bool const & setting ){ fast_ = setting; }

		void set_medium_fast( bool const & setting ){ medium_fast_ = setting; }

		void set_sampler_native_rmsd_screen( bool const & setting ){ sampler_native_rmsd_screen_ = setting; }

		void set_o2star_screen( bool const & setting ){ o2star_screen_ = setting; }

		void set_verbose( bool const & setting ){ verbose_ = setting; }

		void set_distinguish_pucker( bool const & setting ){ distinguish_pucker_ = setting; }

		void set_finer_sampling_at_chain_closure( bool const & setting ){ finer_sampling_at_chain_closure_ = setting; }

		void set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

		void set_allow_syn_pyrimidine( bool const & setting ){ allow_syn_pyrimidine_ = setting; }

		void set_extra_syn_chi_rotamer( bool const & setting ){ extra_syn_chi_rotamer_ = setting; }

		void set_extra_anti_chi_rotamer( bool const & setting ){ extra_anti_chi_rotamer_ = setting; }

		void set_use_phenix_geo( bool const & setting ){ use_phenix_geo_ = setting; }

		void set_centroid_screen( bool const & setting ){ centroid_screen_ = setting; }

		void set_VDW_atr_rep_screen( bool const & setting ){ VDW_atr_rep_screen_ = setting; }

		void set_force_centroid_interaction( bool const & setting ){ force_centroid_interaction_ = setting; }

		void set_choose_random( bool const & setting ){ choose_random_ = setting; }

		void set_skip_sampling( bool const & setting ){ skip_sampling_ = setting; }

		void set_perform_minimize( bool const & setting ){ perform_minimize_ = setting; }

		void set_minimize_and_score_sugar( bool const & setting ){ minimize_and_score_sugar_ = setting; }

		void set_minimize_and_score_native_pose( bool const & setting ){ minimize_and_score_native_pose_ = setting; }

		void set_rm_virt_phosphate( bool const & setting ){ rm_virt_phosphate_ = setting; }

		void set_VDW_rep_screen_info( utility::vector1< std::string > const & setting ){ VDW_rep_screen_info_ = setting; }

		void set_VDW_rep_alignment_RMSD_CUTOFF( core::Real const & setting ){ VDW_rep_alignment_RMSD_CUTOFF_ = setting; }

		void set_output_pdb( bool const & setting ){ output_pdb_ = setting; }

		void set_output_minimized_pose_data_list( bool const & setting ){ output_minimized_pose_data_list_ = setting; }

		void set_fixed_res( utility::vector1< Size > const & setting ){ fixed_res_ = setting; }

	public:
		StepWiseRNA_JobParametersOP
		setup_job_parameters_for_erraser( utility::vector1< Size > moving_res, core::pose::Pose const & pose );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010
		core::pose::PoseCOP native_pose_;

		utility::vector1< core::Size > const moving_res_;
		utility::vector1< core::Size > fixed_res_;
		core::scoring::ScoreFunctionOP scorefxn_;
		std::string silent_file_;
		core::Size sampler_num_pose_kept_;
		core::Size num_pose_minimize_;
		core::Size nstruct_;
		core::Size num_sampled_;
		core::Real sampler_native_screen_rmsd_cutoff_;
		core::Real cluster_rmsd_;
		core::Real native_edensity_score_cutoff_;
		bool fast_;
		bool medium_fast_;
		bool sampler_native_rmsd_screen_;
		bool o2star_screen_;
		bool verbose_;
		bool distinguish_pucker_;
		bool finer_sampling_at_chain_closure_;
		bool PBP_clustering_at_chain_closure_;
		bool allow_syn_pyrimidine_;
		bool extra_syn_chi_rotamer_;
		bool extra_anti_chi_rotamer_;
		bool use_phenix_geo_;
		bool centroid_screen_;
		bool VDW_atr_rep_screen_;
		bool force_centroid_interaction_;
		bool choose_random_;
		bool skip_sampling_;
		bool perform_minimize_;
		bool minimize_and_score_sugar_;
		bool minimize_and_score_native_pose_;
		bool rm_virt_phosphate_;

		utility::vector1< std::string > VDW_rep_screen_info_;
		core::Real VDW_rep_alignment_RMSD_CUTOFF_;
		bool output_pdb_;
		bool output_minimized_pose_data_list_;

	};

} //rna
} //swa
} //protocols

#endif
