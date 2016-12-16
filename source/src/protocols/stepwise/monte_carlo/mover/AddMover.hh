// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddMover.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_monte_carlo_AddMover_hh
#define INCLUDED_protocols_stepwise_monte_carlo_AddMover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_TorsionMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.fwd.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class AddMover: public protocols::moves::Mover {
public:

	AddMover( core::scoring::ScoreFunctionCOP scorefxn = 0 );

	//destructor
	~AddMover();

public:

	using protocols::moves::Mover::apply;

	void
	apply( core::pose::Pose & pose, core::Size const res_to_add_in_full_model_numbering, core::Size const res_to_build_off_in_full_model_numbering );

	void
	apply( core::pose::Pose & pose, StepWiseMove const & swa_move );

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void set_kT( core::Real const & setting ){ kT_ = setting; }

	void set_sample_range_small( core::Real const setting ){ sample_range_small_ = setting; }

	void set_sample_range_large( core::Real const setting ){ sample_range_large_ = setting; }

	void set_internal_cycles( core::Size const setting ){ internal_cycles_ = setting; }

	void set_sample_pH( bool const setting ){ sample_pH_ = setting; }

	void set_presample_added_residue( core::Size const setting ){ presample_added_residue_ = setting; }

	void set_presample_by_swa( core::Size const setting ){ presample_by_swa_ = setting; }

	void set_start_added_residue_in_aform( core::Size const setting ){ start_added_residue_in_aform_ = setting; }

	void set_minimize_single_res( core::Size const setting ){ minimize_single_res_ = setting; }

	void set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler );

private:

	void sample_by_swa( core::pose::Pose & pose, core::Size const res_to_add  ) const;

	void sample_by_monte_carlo_internal( core::pose::Pose & pose ) const;

	void
	do_append( core::pose::Pose & pose );

	void
	do_prepend( core::pose::Pose & pose );

	void
	append_other_pose( core::pose::Pose & pose, core::Size const offset,
		core::Size const other_pose_idx );

	void
	prepend_other_pose( core::pose::Pose & pose, core::Size const offset,
		core::Size const other_pose_idx );

	void
	append_residue( core::pose::Pose & pose, core::Size const offset );

	void
	prepend_residue( core::pose::Pose & pose, core::Size const offset );

	bool
	check_same_chain( core::pose::Pose const & pose, core::Size const res_to_add_in_full_model_numbering, core::Size const res_to_build_off_in_full_model_numbering );

	bool
	check_same_chain( core::pose::Pose const & pose );

	void
	setup_initial_torsions( core::pose::Pose & pose );

	core::Size
	get_add_res( StepWiseMove const & swa_move, core::pose::Pose const & pose ) const;

	void
	setup_initial_jump( core::pose::Pose & pose );

	core::conformation::ResidueOP
	create_residue_to_add( core::pose::Pose const & pose );

private:

	core::chemical::ResidueTypeSetCAP rsd_set_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool presample_added_residue_;
	bool presample_by_swa_;
	bool minimize_single_res_;
	bool start_added_residue_in_aform_;
	core::Size internal_cycles_;
	bool sample_pH_ = false;

	rna::RNA_TorsionMoverOP rna_torsion_mover_;
	core::Real sample_range_small_;
	core::Real sample_range_large_;
	core::Real kT_;

	protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler_;

	core::Size suite_num_, nucleoside_num_;
	core::Size res_to_add_in_full_model_numbering_, res_to_build_off_in_full_model_numbering_;
	StepWiseMove swa_move_;

};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
