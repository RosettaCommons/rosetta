// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Clusterer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das, Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_Clusterer_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_Clusterer_HH

#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

#include <protocols/stepwise/sampling/rna/legacy/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>

#include <map>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {

	class SlicedPoseWorkingParameters: public utility::pointer::ReferenceCount{

		public:

			 SlicedPoseWorkingParameters():
				is_setup_( false )
			{
			}

			virtual ~SlicedPoseWorkingParameters(); // auto-removing definition from header{};

			void
			setup( protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP & working_parameters );

			core::pose::Pose
			create_sliced_pose( core::pose::Pose const & working_pose );

		public:
			utility::vector1< core::Size > sliced_pose_best_alignment; //check
			utility::vector1 < core::Size > sliced_pose_calc_rms_res;
			std::map< core::Size, core::Size > sliced_pose_full_to_sub;
			std::map< core::Size, bool > sliced_pose_is_prepend_map;

 	 	private:
			utility::vector1< bool > is_sliced_res_;
			utility::vector1 < core::Size > working_to_sliced_res_map_;
			utility::vector1 < core::Size > sliced_to_working_res_map_;
			utility::vector1< std::pair < core::Size, core::Size > > delete_res_range_list_;
			bool is_setup_;

	};

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
	class Cluster_Member{

		public:

		Cluster_Member():
		ID( 0 ),
		RMSD( 0.0 ),
		score( 9999.99 )
		{
		}

		~Cluster_Member(){};

		public:

		core::Size ID;
		core::Real RMSD;
		core::Real score;

	};


  /////////////////////////////////////////////////////////////////////////////////////////////////

  class StepWiseRNA_Clusterer: public utility::pointer::ReferenceCount {
  public:

    //constructor!
		StepWiseRNA_Clusterer( utility::vector1< std::string > const & silent_files_in );

		StepWiseRNA_Clusterer( std::string const & silent_file_in );

		StepWiseRNA_Clusterer(  core::io::silent::SilentFileDataOP & sfd );

    //destructor -- necessary?
    virtual ~StepWiseRNA_Clusterer();

    /// @brief Filter a list of poses by score.


		void set_max_decoys( core::Size const & setting ){ max_decoys_ = setting; }

		void set_score_diff_cut( core::Real const & setting ){ score_diff_cut_ = setting; }

		void set_perform_score_diff_cut( bool const & setting ){ perform_score_diff_cut_ = setting; }

		void set_cluster_radius( core::Real const & setting ){ whole_struct_cluster_radius_ = setting; }

		void set_suite_cluster_radius( core::Real const & setting ){ suite_cluster_radius_ = setting; }

		void set_loop_cluster_radius( core::Real const & setting ){ loop_cluster_radius_ = setting; }

		void set_rename_tags( core::Real const & setting ){ rename_tags_ = setting; }

		void set_distinguish_pucker( core::Real const & setting ){ distinguish_pucker_ = setting; }

		void set_add_lead_zero_to_tag( core::Real const & setting ){ add_lead_zero_to_tag_ = setting; }

		void set_working_parameters( protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP & working_parameters );

		void set_working_parameters_exist( bool const working_parameters_exist );

		void set_quick_alignment( bool const & setting ){ quick_alignment_ = setting; }

		void set_align_only_over_base_atoms( bool const & setting ){ align_only_over_base_atoms_ = setting; }

		void set_optimize_memory_usage(  bool const & setting ){ optimize_memory_usage_ = setting; }

		void set_keep_pose_in_memory(  bool const & setting ){ keep_pose_in_memory_ = setting; }

		void set_two_stage_clustering(  bool const & setting ){ two_stage_clustering_ = setting; }

		void set_verbose(  bool const & setting ){ verbose_ = setting; }


		void cluster();

		void
		output_silent_file( std::string const & silent_file );

		void
		recalculate_rmsd_and_output_silent_file( std::string const & silent_file,
				                                    protocols::stepwise::sampling::rna::legacy::StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup,
																					bool const write_score_only );

		void
		get_best_neighboring_shift_RMSD_and_output_silent_file( std::string const & silent_file );

		void
		set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

		void
		set_skip_clustering( bool const & setting ){ skip_clustering_ = setting; }

		void
		set_filter_virtual_res_list( utility::vector1 < core::Size > const & setting ){ filter_virtual_res_list_ = setting; }

		void
		set_perform_VDW_rep_screen( bool const & setting ){ perform_VDW_rep_screen_ = setting; }

		void
		set_perform_filters( bool const & setting ){ perform_filters_ = setting; }

		void
		set_min_num_south_sugar_filter( Size const & setting ){ min_num_south_sugar_filter_ = setting; }

		void
		set_VDW_rep_screen_info( utility::vector1< std::string > const & setting ){ VDW_rep_screen_info_ = setting; }

		void
		set_user_input_VDW_bin_checker( protocols::stepwise::sampling::rna::checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }

		void
		set_full_length_loop_rmsd_clustering( bool const & setting ){ full_length_loop_rmsd_clustering_ = setting; }

		void
		set_ignore_FARFAR_no_auto_bulge_tag( bool const & setting ){ ignore_FARFAR_no_auto_bulge_tag_ = setting; }

		void
		set_ignore_FARFAR_no_auto_bulge_parent_tag( bool const & setting ){ ignore_FARFAR_no_auto_bulge_parent_tag_ = setting; }

		void
		set_ignore_unmatched_virtual_res( bool const & setting ){ ignore_unmatched_virtual_res_ = setting; }

		void
		set_output_pdb( bool const setting ){ output_pdb_ = setting; }

  private:

		void
		initialize_parameters_and_input();

		void
		create_silent_file_and_tag_list();

		void
		do_some_clustering();

		void
		two_stage_clustering();

		void
		create_large_cluster_centers_member_list();


		bool
		is_old_individual_suite_cluster( core::pose::Pose const & current_pose,
                                    core::pose::Pose const & cluster_center_pose,
                                    utility::vector1 < core::Size > const & calc_rms_res,
																	std::map< core::Size, core::Size > const & full_to_sub,
																	std::map< core::Size, bool > const & is_prepend_map,
																	core::Real const & cluster_radius ) const;


		core::pose::PoseOP
		get_poseOP( Size const n );

		void
		setup_fail_triangle_inequailty_list( core::pose::Pose & current_pose, std::string const & tag, utility::vector1< bool > & fail_triangle_inequality_list );


		bool
		is_new_cluster_center_with_working_parameters( core::pose::PoseOP const & pose_op, std::string const & tag );

		bool
		check_for_closeness_without_working_parameters( core::pose::PoseOP const & pose_op );


		bool
		check_for_closeness( core::pose::PoseOP const & pose_op, std::string const & tag );

		utility::vector1< core::Size > const &
		get_act_alignment_res()  const;

		utility::vector1 < core::Size > const &
		get_act_calc_rms_res()	 const;

		std::map< core::Size, core::Size > const &
		get_act_full_to_sub()	 const;

		std::map< core::Size, bool > const &
		get_act_is_prepend_map() const;

		void
		initialize_quick_alignment_pose();

		void
		initialize_max_memory_pose_num();

		void
		align_to_quick_alignment_pose( core::pose::Pose & pose, std::string const & tag ) const;

		void
		initialize_VDW_rep_checker();

		void
		create_tags_map();

		bool
		pass_FARFAR_no_auto_bulge_filter( core::io::silent::SilentStructOP const & silent_struct ) const;


  private:

		utility::vector1< std::string > silent_files_;
		core::import_pose::pose_stream::SilentFilePoseInputStreamOP input_;

		utility::vector1< core::pose::PoseOP > pose_output_list_;
		utility::vector1< std::string > tag_output_list_;
		utility::vector1< core::io::silent::SilentStructOP > silent_struct_output_list_;

		Size max_decoys_;
		core::Real score_diff_cut_;
		bool perform_score_diff_cut_;

		core::Real whole_struct_cluster_radius_;
		core::Real suite_cluster_radius_;
		core::Real loop_cluster_radius_;

		bool rename_tags_;
		protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP working_parameters_;
		bool working_parameters_exist_;
		bool distinguish_pucker_;
		bool add_lead_zero_to_tag_;
		bool quick_alignment_;
		bool align_only_over_base_atoms_;
		bool optimize_memory_usage_;
		SlicedPoseWorkingParameters sliced_pose_job_params_;

		bool verbose_;
		bool keep_pose_in_memory_;
		bool keep_pose_in_memory_hydrid_;
		Size max_memory_pose_num_;

		bool two_stage_clustering_;
		bool use_triangle_inequality_;
		core::pose::Pose first_pose_;
		utility::vector1< utility::vector1< Cluster_Member > > cluster_centers_neighbor_list_;
		utility::vector1< core::pose::PoseOP > large_cluster_pose_list_;
		utility::vector1< core::Size > all_pose_to_output_pose_ID_map_;
		bool PBP_clustering_at_chain_closure_;

		bool quick_alignment_pose_is_intialized_;
		core::pose::Pose quick_alignment_pose_;
		std::string quick_alignment_tag_;
		bool skip_clustering_;

		protocols::stepwise::sampling::rna::checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;
		bool perform_VDW_rep_screen_;
		bool perform_filters_;
		Size min_num_south_sugar_filter_;
		utility::vector1< std::string > VDW_rep_screen_info_;
		bool full_length_loop_rmsd_clustering_;

		utility::vector1 < core::Size > filter_virtual_res_list_;

		std::map< std::string, bool > current_tags_map_;
		std::map< std::string, bool > parent_tags_map_;

		bool ignore_FARFAR_no_auto_bulge_parent_tag_;
		bool ignore_FARFAR_no_auto_bulge_tag_;
		bool ignore_unmatched_virtual_res_;
		bool output_pdb_;
		core::chemical::ResidueTypeSetCAP rsd_set_;

  };

} //rna
} //sampling
} //stepwise
} //protocols

#endif
