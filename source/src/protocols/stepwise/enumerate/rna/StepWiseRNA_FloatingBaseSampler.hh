// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/StepWiseRNA_FloatingBaseSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSampler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSampler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/enumerate/rna/FloatingBaseClasses.hh>
#include <protocols/stepwise/enumerate/rna/SugarModeling.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_PoseSelection.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/AtrRepScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/ChainBreakScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/ChainClosableScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/O2PrimePacker.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {

	class StepWiseRNA_FloatingBaseSampler: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseRNA_FloatingBaseSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_FloatingBaseSampler();

		virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

	public:

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_scorefxn( scoring::ScoreFunctionOP const & scorefxn );

		utility::vector1< pose::PoseOP > &		get_pose_list();

		void
		set_base_centroid_screener( screener::StepWiseRNA_BaseCentroidScreenerOP & screener );

		void
		set_user_input_VDW_bin_screener( screener::StepWiseRNA_VDW_BinScreenerOP const & user_input_VDW_bin_screener );

		void
		set_anchor_sugar_modeling( SugarModeling const & anchor_sugar_modeling );

		void
		set_options( StepWiseRNA_ModelerOptionsCOP options );

	private:

		void
		initialize_poses_and_stubs_and_screeners( pose::Pose & pose  );

		void
		reinstantiate_backbone_and_add_constraint_at_moving_res(	pose::Pose & pose, Size const & five_prime_chain_break_res );

		void
		initialize_euler_angle_grid_parameters();

		void
		initialize_xyz_grid_parameters();

		void
		initialize_rigid_body_sampler( pose::Pose const & pose );

		bool
		break_early_for_integration_tests();

		bool
		screen_anchor_sugar_conformation( pose::Pose & pose, std::string & tag );

		kinematics::Stub
		get_reference_stub( pose::Pose const & pose ) const;

		void
		initialize_other_residues_base_list( pose::Pose const & pose );

		void
		update_base_bin_map( Base_bin const & base_bin );

		void
		update_base_bin_map( utility::vector1< Real > const & rigid_body_values );

		void
		output_count_data();

		void
		instantiate_moving_sugar_and_o2prime( pose::Pose & pose );

		void
		virtualize_moving_sugar_and_o2prime( pose::Pose & pose );

		bool
		check_moving_sugar( pose::Pose & pose );

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

		protocols::rotamer_sampler::rigid_body::RigidBodyRotamerOP sampler_;

		StepWiseRNA_CountStruct count_data_;
		std::string silent_file_;
		bool native_rmsd_screen_;
		Size num_pose_kept_to_use_;

		StepWiseRNA_ModelerOptionsCOP options_;

		SugarModeling anchor_sugar_modeling_;

		screener::AtrRepScreenerOP atr_rep_screener_;
		screener::AtrRepScreenerOP atr_rep_screener_with_instantiated_sugar_;
		utility::vector1< screener::AtrRepScreenerOP > atr_rep_screeners_for_anchor_sugar_models_;
		screener::StepWiseRNA_VDW_BinScreenerOP VDW_bin_screener_;
		screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;
		screener::ChainClosableScreenerOP chain_closable_to_distal_screener_, chain_closable_to_anchor_screener_;
		screener::ChainBreakScreenerOP chain_break_to_distal_screener_, chain_break_to_anchor_screener_;
		screener::StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;

		pose::PoseOP screening_pose_, sugar_screening_pose_;
		utility::vector1 < conformation::ResidueOP > moving_rsd_at_origin_list_, screening_moving_rsd_at_origin_list_, sugar_screening_moving_rsd_at_origin_list_;
		StepWiseRNA_PoseSelectionOP pose_selection_;

		O2PrimePackerOP o2prime_packer_;

		kinematics::Stub reference_stub_;
		utility::vector1 < core::kinematics::Stub > other_residues_base_list_;
		std::map< Base_bin, int, compare_base_bin > base_bin_map_;

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
	};

} //rna
} //enumerate
} //stepwise
} //protocols

#endif
