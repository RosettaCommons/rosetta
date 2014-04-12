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
#include <protocols/stepwise/sampling/protein/StepWiseProteinJobParameters.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>


using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinCCD_Closer: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinCCD_Closer(rotamer_sampler::RotamerSizedOP sampler,
															protocols::stepwise::sampling::protein::StepWiseProteinJobParametersCOP job_parameters );

    //destructor -- necessary?
    ~StepWiseProteinCCD_Closer();

    /// @brief Apply the minimizer to one pose
    virtual void apply( pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		utility::vector1< utility::vector1< Real > > const & main_chain_torsion_set_lists() const;
		utility::vector1< id::TorsionID > const & which_torsions() const;

		void
		set_ccd_close_res( Size const value ){ ccd_close_res_ = value;}

		void
		set_working_moving_res_list( utility::vector1< Size > const & setting ){ moving_residues_ = setting; }

		void set_choose_random( bool const & setting ){ choose_random_ = setting; }
		bool choose_random() const{ return choose_random_; }

		void set_num_random_samples( Size const & setting ){ num_random_samples_ = setting; }
		Size num_random_samples() const{ return num_random_samples_; }

	private:

		bool
		CCD_loop_close( pose::Pose & pose );

		void
		CCD_loop_close_sample_omega_recursively( pose::Pose & pose, int const offset );

		void
		setup_torsions( pose::Pose const & pose );

		void
		figure_out_loop( pose::Pose const & pose );

		void
		grab_main_chain_torsion_set_list( pose::Pose const & pose );

		void
		save_phi_psi_omega_over_loop_residues( pose::Pose const & pose );

		void
		restore_phi_psi_omega_over_loop_residues( pose::Pose & pose );

		void
		restore_phi_psi_over_loop_residues( pose::Pose & pose );

		void
		fix_jump_atoms_at_loop_boundaries( pose::Pose & pose );

	private:

		rotamer_sampler::RotamerSizedOP sampler_;

		utility::vector1< Size > working_bridge_res_;
		utility::vector1< Size > moving_residues_;
		utility::vector1< bool > is_pre_proline_;
		Size ccd_close_res_;

		loops::Loop loop_;
		kinematics::MoveMap mm_;

		utility::vector1< id::TorsionID >  which_torsions_;
		utility::vector1< utility::vector1< Real > > main_chain_torsion_sets_for_moving_residues_;

		utility::vector1< Real > main_chain_torsion_set_for_moving_residues_save_;

		bool const verbose_;

		Size pose_count_;

		bool choose_random_;
		Size num_random_samples_;
		Size max_ntries_;

  };

} //protein
} //sampling
} //stepwise
} //protocols

#endif
