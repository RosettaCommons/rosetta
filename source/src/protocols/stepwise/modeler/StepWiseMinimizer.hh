// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/StepWiseMinimizer.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_StepWiseMinimizer_HH
#define INCLUDED_protocols_stepwise_modeler_StepWiseMinimizer_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/modeler/StepWiseMinimizer.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh>
#include <protocols/toolbox/AllowInsert.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/CartesianMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>

namespace protocols {
namespace stepwise {
namespace modeler {

	class StepWiseMinimizer: public moves::Mover {

	public:

		//constructor
		StepWiseMinimizer( utility::vector1< core::pose::PoseOP > const & pose_list,
											 working_parameters::StepWiseWorkingParametersCOP working_parameters,
											 options::StepWiseModelerOptionsCOP options,
											 core::scoring::ScoreFunctionCOP scorefxn);


		//destructor
		~StepWiseMinimizer();

	public:

		virtual void apply( core::pose::Pose & pose );

		virtual std::string get_name() const{ return "StepWiseMinimizer"; }

	public:

		void
		set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){ scorefxn_ = scorefxn; }

    void set_working_pack_res( utility::vector1< core::Size > const & setting ){ working_pack_res_ = setting; }

		utility::vector1< core::pose::PoseOP > pose_list() const { return pose_list_; }

  private:

		void
		initialize_parameters();

		void
		let_neighboring_side_chains_minimize( core::kinematics::MoveMap & mm,
																				  core::pose::Pose const & pose );
		void setup_minimizers();

		void setup_scorefxns( core::pose::Pose const & pose );

		void do_minimize( core::pose::Pose & pose, core::kinematics::MoveMap & mm );

		void
		do_full_minimizing( core::pose::Pose & pose );

		void
		do_clustering( core::pose::Pose & pose );

		void get_move_map_and_allow_insert( core::kinematics::MoveMap & mm, core::pose::Pose const & pose );

		void close_chainbreaks( core::pose::Pose & pose, core::kinematics::MoveMap & mm );

		void setup_all_moving_res( core::pose::Pose const & pose );

		void
		move_side_chain( core::kinematics::MoveMap & mm,
										 core::pose::Pose const & pose,
										 core::Size const j );

		void
		setup_vary_bond_geometry( core::pose::Pose & pose, core::kinematics::MoveMap & mm );

	private:

		utility::vector1< core::pose::PoseOP > pose_list_;
		options::StepWiseModelerOptionsCOP options_;
		core::scoring::ScoreFunctionCOP scorefxn_;
		core::scoring::ScoreFunctionOP minimize_scorefxn_, final_scorefxn_;

		Size num_pose_minimize_;
    utility::vector1< core::Size > const working_moving_res_;
    utility::vector1< core::Size > working_fixed_res_;
    utility::vector1< core::Size > working_calc_rms_res_;
    utility::vector1< core::Size > working_pack_res_;

		core::optimization::CartesianMinimizerOP cartesian_minimizer_;
		core::optimization::AtomTreeMinimizerOP atom_tree_minimizer_;
		core::optimization::MinimizerOptionsOP minimizer_options_;
		bool const allow_virtual_o2prime_hydrogens_;

		protein::loop_close::StepWiseProteinCCD_CloserOP protein_ccd_closer_;

		protocols::toolbox::AllowInsertOP allow_insert_; // a atom-centric version of the DOF-centric movemap.
	};

} //modeler
} //stepwise
} //protocols

#endif
