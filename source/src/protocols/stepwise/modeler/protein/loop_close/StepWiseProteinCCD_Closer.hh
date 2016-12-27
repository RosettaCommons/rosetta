// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SWA_ProteinCCD_Closer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinCCD_Closer_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinCCD_Closer_HH

#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.fwd.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace loop_close {

/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseProteinCCD_Closer: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseProteinCCD_Closer( protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters );

	//destructor -- necessary?
	virtual ~StepWiseProteinCCD_Closer();

	void
	init( core::pose::Pose & pose );

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	void get_closure_solution( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	utility::vector1< core::id::TorsionID > const & which_torsions() const;

	void
	set_ccd_close_res( core::Size const value ){ ccd_close_res_ = value;}

	void
	set_working_moving_res_list( utility::vector1< Size > const & setting ){ moving_residues_ = setting; }

	utility::vector1< core::Real >
	grab_main_chain_torsion_set_list( core::pose::Pose const & pose );

	utility::vector1< core::Real > const &
	main_chain_torsion_set() const;

	utility::vector1< core::Real > const &
	main_chain_torsion_set_save() const;

	core::Size ntries() const { return ntries_; }

	bool closed_loop() const { return closed_loop_; }

private:

	bool
	CCD_loop_close( core::pose::Pose & pose );

	void
	CCD_loop_close_sample_omega_recursively( core::pose::Pose & pose, int const offset );

	void
	setup_torsions();

	void
	figure_out_loop( core::pose::Pose const & pose );

	void
	figure_out_movemap();

	void
	save_phi_psi_omega_over_loop_residues( core::pose::Pose const & pose );

	void
	restore_phi_psi_omega_over_loop_residues( core::pose::Pose & pose );

	void
	restore_phi_psi_over_loop_residues( core::pose::Pose & pose );

	void
	fix_jump_atoms_at_loop_boundaries( core::pose::Pose & pose );

	Size
	check_for_unique_cutpoint_flanked_by_bridge_res( core::pose::Pose const & pose );

	Size
	check_for_unique_cutpoint( core::pose::Pose const & pose );

private:

	loops::loop_closure::ccd::CCDLoopClosureMoverOP ccd_loop_closure_mover_;

	utility::vector1< Size > working_bridge_res_;
	utility::vector1< Size > moving_residues_;
	utility::vector1< bool > is_pre_proline_;
	core::Size ccd_close_res_;

	loops::Loop loop_;
	core::kinematics::MoveMapOP mm_;

	utility::vector1< core::id::TorsionID >  which_torsions_;

	utility::vector1< core::Real > main_chain_torsion_set_;
	utility::vector1< core::Real > main_chain_torsion_set_save_;

	bool closed_loop_;
	Size ntries_;
};

} //loop_close
} //protein
} //modeler
} //stepwise
} //protocols

#endif
