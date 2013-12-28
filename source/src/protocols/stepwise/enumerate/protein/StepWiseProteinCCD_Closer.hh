// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ProteinCCD_Closer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinCCD_Closer_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinCCD_Closer_HH

#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/stepwise/enumerate/protein/StepWiseJobParameters.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.fwd.hh>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinCCD_Closer: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinCCD_Closer(
															 protocols::stepwise::enumerate::protein::sample_generators::StepWisePoseSampleGeneratorOP sample_generator,
															 protocols::stepwise::enumerate::protein::StepWiseJobParametersOP job_parameters );

    //destructor -- necessary?
    ~StepWiseProteinCCD_Closer();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;


		utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists() const;
		utility::vector1< core::id::TorsionID > const & which_torsions() const;

		void
		set_ccd_close_res( Size const value ){ ccd_close_res_ = value;}

	private:

		bool
		CCD_loop_close( core::pose::Pose & pose );

		void
		CCD_loop_close_sample_omega_recursively( core::pose::Pose & pose, int const offset );

		void
		setup_torsions( core::pose::Pose const & pose );

		void
		figure_out_loop( core::pose::Pose const & pose );

		void
		grab_main_chain_torsion_set_list( core::pose::Pose const & pose );

		void
		save_phi_psi_omega_over_loop_residues( core::pose::Pose const & pose );

		void
		restore_phi_psi_omega_over_loop_residues( core::pose::Pose & pose );

		void
		restore_phi_psi_over_loop_residues( core::pose::Pose & pose );

	private:

		protocols::stepwise::enumerate::protein::sample_generators::StepWisePoseSampleGeneratorOP sample_generator_;

		utility::vector1< Size > working_bridge_res_;
		utility::vector1< Size > moving_residues_;
		utility::vector1< bool > is_pre_proline_;
		Size ccd_close_res_;

		loops::Loop loop_;
		core::kinematics::MoveMap mm_;

		utility::vector1< core::id::TorsionID >  which_torsions_;
		utility::vector1< utility::vector1< core::Real > > main_chain_torsion_sets_for_moving_residues_;

		utility::vector1< core::Real > main_chain_torsion_set_for_moving_residues_save_;

		bool const verbose_;

		Size pose_count_;
  };

} //protein
} //enumerate
} //stepwise
} //protocols

#endif
