// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
/// @brief
/// @details
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_VirtualSugarSampler_HH
#define INCLUDED_protocols_stepwise_rna_VirtualSugarSampler_HH

#include <protocols/moves/MoverForPoseList.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009

/*
Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using

using namespace core;
using namespace core::pose;
*/

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {


class VirtualSugarSampler: public protocols::moves::MoverForPoseList {

public:

	//constructor
	VirtualSugarSampler( working_parameters::StepWiseWorkingParametersCOP & working_parameters, SugarModeling & sugar_modeling );

	//destructor
	~VirtualSugarSampler();

	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual void
	apply( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	void set_tag( std::string const & setting ) { tag_ = setting; }

	void set_use_phenix_geo( bool const & setting ) { use_phenix_geo_ = setting; }

	void set_legacy_mode( bool const & setting ) { legacy_mode_ = setting; }

	void set_choose_random( bool const & setting ) { choose_random_ = setting; }

	void set_keep_base_fixed( bool const & setting ) { keep_base_fixed_ = setting; }

	void set_do_minimize( bool const & setting ) { do_minimize_ = setting; }

	void set_do_screens( bool const & setting ) { do_screens_ = setting; }

	void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

	void set_virtual_sugar_is_from_prior_step( bool const & setting ) { virtual_sugar_is_from_prior_step_ = setting; }

	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

private:

	void
	setup_sugar_conformations( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & pose );

	void
	minimize_sugar( core::pose::Pose & pose_with_sugar );

	void
	get_sugar_setup_scorefxns( core::scoring::ScoreFunctionOP & sugar_scorefxn, core::scoring::ScoreFunctionOP & sugar_scorefxn_without_ch_bond, core::scoring::ScoreFunctionOP & rescaled_sugar_score_fxn_without_ch_bond ) const;

	void
	do_chain_closure_modeler( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & viewer_pose );

	void
	bulge_chain_closure( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & viewer_pose );

	void
	bulge_chain_closure_complete( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & viewer_pose );

	void
	bulge_chain_closure_legacy( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & viewer_pose );

	void
	bulge_chain_minimize_legacy( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose & viewer_pose );

	void
	reinstantiate_backbone_at_moving_res( core::pose::Pose & pose, core::Size const rebuild_res,
		core::Size const five_prime_chain_break_res );

	void
	initialize_pose_variants_for_chain_closure( utility::vector1< core::pose::PoseOP > & pose_list );

	void
	restore_pose_variants_after_chain_closure( utility::vector1< core::pose::PoseOP > & pose_list );

	bool
	fast_full_atom_VDW_repulsion_screen( core::pose::Pose const & pose, core::Size const res_1, core::Size const res_2, bool const is_prepend );

	void
	setup_VDW_bin_checker( core::pose::Pose const & input_pose );

	void
	virtualize_distal_partition( core::pose::Pose & input_pose );

	void
	reinstantiate_distal_partition( utility::vector1< core::pose::PoseOP > & final_pose_list );

	void
	reinstantiate_distal_partition( core::pose::Pose & current_pose );

	void
	reinstate_original_constraints( utility::vector1< core::pose::PoseOP >  & pose_list );

private:

	working_parameters::StepWiseWorkingParametersCOP working_parameters_;
	SugarModeling & sugar_modeling_; // trick -- inputs some modeling info, and holds poses_list as output.
	std::string tag_;
	bool use_phenix_geo_;
	bool keep_base_fixed_;
	bool choose_random_;
	bool do_minimize_;
	bool do_screens_;
	bool integration_test_mode_;
	bool virtual_sugar_is_from_prior_step_;
	bool legacy_mode_;
	bool const do_chain_closure_;
	bool const first_minimize_with_fixed_base_;
	Size const max_tries_for_random_sugar_setup_;
	bool sugar_setup_success_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	checker::RNA_VDW_BinCheckerOP VDW_bin_checker_;

	utility::vector1 < core::Size > distal_partition_res_;
	utility::vector1 < core::Size > already_virtualized_res_list_;
	bool moving_phosphate_virtualized_;
	core::pose::PoseOP pose_with_original_terminal_phosphates_;

	core::scoring::constraints::ConstraintSetOP original_constraint_set_;

};


} //sugar
} //rna
} //modeler
} //stepwise
} //protocols

#endif
