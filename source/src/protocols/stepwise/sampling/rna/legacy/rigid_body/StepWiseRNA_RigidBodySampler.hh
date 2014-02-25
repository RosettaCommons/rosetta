// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/rigid_body/StepWiseRNA_RigidBodySampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_RigidBodySampler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_RigidBodySampler_HH

#include <protocols/moves/MoverForPoseList.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/rigid_body/FloatingBaseClasses.hh>
#include <protocols/stepwise/sampling/rna/sugar/SugarModeling.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosureChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {
namespace rigid_body {

	class StepWiseRNA_RigidBodySampler: public protocols::moves::MoverForPoseList {

	public:

		//constructor
		StepWiseRNA_RigidBodySampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_RigidBodySampler();

		virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		using MoverForPoseList::apply;

	public:

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_scorefxn( scoring::ScoreFunctionOP const & scorefxn );

		utility::vector1< pose::PoseOP > &		get_pose_list();

		void
		set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker );

		void
		set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker );

		void
		set_anchor_sugar_modeling( sugar::SugarModeling const & anchor_sugar_modeling );

		void
		set_options( StepWiseRNA_ModelerOptionsCOP options );

	private:

		void
		initialize_poses_and_stubs_and_checkers( pose::Pose & pose  );

		void
		reinstantiate_backbone_and_add_constraint_at_moving_res(	pose::Pose & pose, Size const & five_prime_chain_break_res );

		void
		initialize_euler_angle_grid_parameters();

		void
		initialize_xyz_grid_parameters();

		void
		initialize_rigid_body_sampler( pose::Pose const & pose );

		rotamer_sampler::rigid_body::RigidBodyRotamerOP
		define_rigid_body_rotamer( pose::Pose const & pose );

		void
		initialize_screeners( pose::Pose & pose );

		void
		initialize_residue_level_screeners( pose::Pose & pose );

		void
		initialize_pose_level_screeners( pose::Pose & pose );


		kinematics::Stub
		get_reference_stub( pose::Pose const & pose ) const;

		void
		update_base_bin_map( BaseBin const & base_bin );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010
		Size const moving_res_; // Might not corresponds to user input.
		Size const moving_suite_; // dofs betweeen this value and value+1 actually move.
		bool const is_prepend_;
		bool const is_internal_; // no cutpoints before or after moving_res.
		Size const gap_size_; /* If this is zero or one, need to screen or closable chain break */
		Size const gap_size_to_anchor_;
		Size const five_prime_chain_break_res_;
		Size const chain_break_reference_res_;
		Size const reference_res_; //the last static_residues that this attach to the moving residues
		Size const floating_base_five_prime_chain_break_;
		Size const floating_base_three_prime_chain_break_;
		bool const is_dinucleotide_;
		bool const close_chain_to_distal_;
		bool const close_chain_to_anchor_;

		utility::vector1< pose::PoseOP > pose_list_;

		scoring::ScoreFunctionOP scorefxn_;

		protocols::rotamer_sampler::rigid_body::RigidBodyRotamerWithResidueAlternativesOP sampler_;
		//		protocols::rotamer_sampler::rigid_body::RigidBodyRotamerWithResidueListOP sampler_;
		utility::vector1< screener::StepWiseScreenerOP > screeners_;
		screener::TagDefinitionOP tag_definition_;

		StepWiseRNA_CountStruct count_data_;
		std::string silent_file_;
		bool native_rmsd_screen_;
		Size num_pose_kept_to_use_;

		StepWiseRNA_ModelerOptionsCOP options_;

		sugar::SugarModeling anchor_sugar_modeling_;

		checker::AtrRepCheckerOP atr_rep_checker_;
		checker::AtrRepCheckerOP atr_rep_checker_with_instantiated_sugar_;
		utility::vector1< checker::AtrRepCheckerOP > atr_rep_checkers_for_anchor_sugar_models_;
		checker::RNA_VDW_BinCheckerOP VDW_bin_checker_;
		checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;
		checker::ChainClosableGeometryCheckerOP chain_closable_geometry_to_distal_checker_, chain_closable_geometry_to_anchor_checker_;
		checker::ChainClosureCheckerOP chain_break_to_distal_checker_, chain_break_to_anchor_checker_;
		checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;

		pose::PoseOP screening_pose_, sugar_screening_pose_;
		utility::vector1< pose::PoseOP > pose_at_origin_list_;
		utility::vector1 < conformation::ResidueOP > moving_rsd_at_origin_list_, screening_moving_rsd_at_origin_list_, sugar_screening_moving_rsd_at_origin_list_;

		StepWiseRNA_PoseSelectionOP pose_selection_;

		o2prime::O2PrimePackerOP o2prime_packer_;
		phosphate::MultiPhosphateSamplerOP phosphate_sampler_;

		kinematics::Stub reference_stub_, moving_res_base_stub_;
		core::conformation::ResidueOP moving_rsd_at_origin_;
		utility::vector1 < core::conformation::ResidueOP > moving_rsd_at_origin_from_sampler_;
		BaseBinMap base_bin_map_;

		Size max_ntries_; // for choose_random

		int euler_angle_bin_min_;
		int euler_angle_bin_max_;
		Real euler_angle_bin_size_;
		int euler_z_bin_min_;
		int euler_z_bin_max_;
		Real euler_z_bin_size_;
		int centroid_bin_min_;
		int centroid_bin_max_;
		Real centroid_bin_size_;
		Distance max_distance_;
		Real max_distance_squared_;

		bool try_sugar_instantiation_;
		Distance o2prime_instantiation_distance_cutoff_;
		std::string const extra_tag_;
		bool const rigid_body_sampling_;
		bool const anchor_sugar_screener_legacy_;
		bool const legacy_mode_;
		bool const residue_level_screening_;
		bool const full_pose_level_screening_;
	};

} //rigid_body
} //legacy
} //rna
} //sampling
} //stepwise
} //protocols

#endif
