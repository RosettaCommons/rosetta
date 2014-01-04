// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/StepWiseRNA_Modeler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_Modeler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_Modeler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Modeler.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Minimizer.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/stepwise/enumerate/rna/screener/StepWiseRNA_VDW_BinScreener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>
#include <string>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {

	class StepWiseRNA_Modeler: public protocols::moves::Mover  {

	public:

		//constructor
		StepWiseRNA_Modeler( scoring::ScoreFunctionOP scorefxn );

		//constructor
		StepWiseRNA_Modeler( Size const sample_res, scoring::ScoreFunctionOP scorefxn );

		//destructor
		~StepWiseRNA_Modeler();

		StepWiseRNA_Modeler( StepWiseRNA_Modeler const & src );

		StepWiseRNA_ModelerOP clone_modeler() const;

		StepWiseRNA_Modeler & operator=( StepWiseRNA_Modeler const & src );

	public:

		void set_moving_res_and_reset( Size const moving_res );

		virtual void apply( pose::Pose & pose );

		virtual std::string get_name() const;

		void set_job_parameters( StepWiseRNA_JobParametersCOP job_parameters );

		void set_native_pose( pose::PoseCOP );

		pose::PoseCOP get_native_pose();

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		Size get_num_sampled(){ return num_sampled_; }

		void set_fixed_res( utility::vector1< Size > const & setting ){ fixed_res_ = setting; minimize_res_.clear(); }

		void set_minimize_res( utility::vector1< Size > const & setting ){ minimize_res_ = setting; fixed_res_.clear(); }

		void set_terminal_res( utility::vector1< Size > const & setting ){ terminal_res_ = setting; }

		void set_rmsd_res_list( utility::vector1< Size > const & setting ){ rmsd_res_list_ = setting; }

		void set_syn_chi_res_list( utility::vector1< core::Size > const & setting ){ syn_chi_res_list_ = setting;}

		void set_skip_sampling( bool const & setting );

		void
		set_minimizer_extra_minimize_res( utility::vector1< core::Size > const & setting ){ minimizer_extra_minimize_res_ = setting; }

		StepWiseRNA_JobParametersOP
		setup_job_parameters_for_stepwise_with_full_model_info( utility::vector1< Size > moving_res,
																											 core::pose::Pose & pose );

		StepWiseRNA_JobParametersOP
		setup_job_parameters_for_swa( utility::vector1< Size > moving_res, pose::Pose const & pose );

		void
		output_pose(pose::Pose & pose, std::string const & out_tag, std::string const out_silent_file ) const;

		StepWiseRNA_ModelerOptionsCOP options();

		void
		set_options( StepWiseRNA_ModelerOptionsCOP options );

	private:

		void initialize_variables();

		void initialize_job_parameters_and_root( pose::Pose & pose );

		void
		do_residue_sampling( pose::Pose & pose, utility::vector1< PoseOP > & pose_list );

		bool
		sampling_successful( utility::vector1< PoseOP > & pose_list );

		void
		do_minimizing( pose::Pose & pose, utility::vector1< PoseOP > & pose_list );

		void
		add_to_pose_list( utility::vector1< pose::PoseOP > & pose_list, pose::Pose const & pose, std::string const pose_tag ) const;

		void
		setup_root_based_on_full_model_info( pose::Pose & pose, StepWiseRNA_JobParametersCOP & job_parameters );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_;
		pose::PoseCOP native_pose_;

		scoring::ScoreFunctionOP scorefxn_;
		std::string silent_file_;

		// maybe following should go into options?
		utility::vector1< Size > moving_res_list_;
		utility::vector1< Size > fixed_res_;
		utility::vector1< Size > minimize_res_;
		utility::vector1< Size > rmsd_res_list_;
		utility::vector1< Size > terminal_res_;
		utility::vector1< core::Size > minimizer_extra_minimize_res_;
		utility::vector1< core::Size > syn_chi_res_list_;

		StepWiseRNA_ModelerOptionsCOP options_;

		Size num_sampled_;
		StepWiseRNA_MinimizerOP stepwise_rna_minimizer_;
		kinematics::MoveMapOP minimize_move_map_;
		screener::StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;
		screener::StepWiseRNA_VDW_BinScreenerOP user_input_VDW_bin_screener_;
		toolbox::AllowInsertOP allow_insert_;

		utility::vector1< std::string > VDW_rep_screen_info_;
	};

} //rna
} //enumerate
} //stepwise
} //protocols

#endif
