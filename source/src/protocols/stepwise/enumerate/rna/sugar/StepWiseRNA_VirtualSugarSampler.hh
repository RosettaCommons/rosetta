// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_VirtualSugarSampler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_VirtualSugarSampler_HH

#include <protocols/moves/MoverForPoseList.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/enumerate/rna/sugar/SugarModeling.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009

using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace sugar {


	class StepWiseRNA_VirtualSugarSampler: public protocols::moves::MoverForPoseList {

	public:

		//constructor
		StepWiseRNA_VirtualSugarSampler( StepWiseRNA_JobParametersCOP & job_parameters, SugarModeling & sugar_modeling	);

		//destructor
		~StepWiseRNA_VirtualSugarSampler();

		virtual void apply( Pose & pose_to_visualize );

		virtual void
		apply( utility::vector1< pose::PoseOP > & pose_list,
					 Pose & pose_to_visualize );

		virtual std::string get_name() const;

		void set_tag( std::string const & setting ) { tag_ = setting;	}

		void set_use_phenix_geo( bool const & setting ) { use_phenix_geo_ = setting;	}

		void set_legacy_mode( bool const & setting ) { legacy_mode_ = setting;	}

		void set_choose_random( bool const & setting ) { choose_random_ = setting;	}

		void set_keep_base_fixed( bool const & setting ) { keep_base_fixed_ = setting;	}

		void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

		void set_virtual_sugar_is_from_prior_step( bool const & setting ) { virtual_sugar_is_from_prior_step_ = setting;	}

		void set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	private:

		void
		setup_sugar_conformations( utility::vector1< PoseOP > & pose_list, pose::Pose & pose );

		void
		minimize_sugar( pose::Pose & pose_with_sugar );

		void
		get_sugar_setup_scorefxns( scoring::ScoreFunctionOP & sugar_scorefxn, scoring::ScoreFunctionOP & sugar_scorefxn_without_ch_bond, scoring::ScoreFunctionOP & rescaled_sugar_score_fxn_without_ch_bond ) const;

		void
		do_chain_closure_sampling( utility::vector1< PoseOP > & pose_list,
															 pose::Pose & viewer_pose );

		void
		initialize_pose_variants_for_chain_closure( utility::vector1< pose::PoseOP > & pose_list );

		void
		restore_pose_variants_after_chain_closure( utility::vector1< pose::PoseOP > & pose_list );

		void
		bulge_chain_closure( utility::vector1< PoseOP > & pose_list,
																 pose::Pose & viewer_pose );

		void
		bulge_chain_closure_complete( utility::vector1< PoseOP > & pose_list,
																					pose::Pose & viewer_pose );

		void
		bulge_chain_closure_legacy( utility::vector1< PoseOP > & pose_list,
																				pose::Pose & viewer_pose );

		void
		bulge_chain_minimize_legacy( utility::vector1< PoseOP > & pose_list,
																 Pose & viewer_pose );

		bool
		fast_full_atom_VDW_repulsion_screen( core::pose::Pose const & pose, core::Size const res_1, core::Size const res_2, bool const is_prepend );

		void
		setup_VDW_bin_screener( pose::Pose const & input_pose );

		void
		virtualize_distal_partition( pose::Pose & input_pose );

		void
		reinstantiate_distal_partition( utility::vector1< PoseOP > & final_pose_list );

		void
		reinstantiate_distal_partition( pose::Pose & current_pose );

		void
		reinstate_original_constraints( utility::vector1< pose::PoseOP >  & pose_list );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_;
		SugarModeling & sugar_modeling_; // trick -- inputs some modeling info, and holds pose_list as output.
		std::string tag_;
		bool use_phenix_geo_;
		bool keep_base_fixed_;
		bool choose_random_;
		bool integration_test_mode_;
		bool virtual_sugar_is_from_prior_step_;
		bool legacy_mode_;
		bool const do_chain_closure_;
		bool const first_minimize_with_fixed_base_;
		Size const max_tries_for_random_overall_;
		Size const max_tries_for_random_sugar_setup_;
		bool sugar_setup_success_;
		scoring::ScoreFunctionOP scorefxn_;
		screener::StepWiseRNA_VDW_BinScreenerOP VDW_bin_screener_;

		utility::vector1 < core::Size > distal_partition_pos_;
		utility::vector1 < core::Size > already_virtualized_res_list_;
		bool moving_phosphate_virtualized_;

		core::scoring::constraints::ConstraintSetOP original_constraint_set_;

	};


} //sugar
} //rna
} //enumerate
} //stepwise
} //protocols

#endif
