// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_MultiPoseCloser.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_loop_close_StepWiseProteinCCD_MultiPoseCloser_HH
#define INCLUDED_protocols_stepwise_sampling_protein_loop_close_StepWiseProteinCCD_MultiPoseCloser_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_MultiPoseCloser.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace loop_close {

	class StepWiseProteinCCD_MultiPoseCloser: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseProteinCCD_MultiPoseCloser( working_parameters::StepWiseWorkingParametersCOP working_parameters,
																			  rotamer_sampler::RotamerSizedOP sampler );

		//destructor
		~StepWiseProteinCCD_MultiPoseCloser();

	public:

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const { return "StepWiseProteinCCD_MultiPoseCloser"; }

		utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_sets() const;

		void
		set_ccd_close_res( core::Size const value ){ ccd_close_res_ = value;}

		void
		set_working_moving_res_list( utility::vector1< core::Size > const & setting ){ moving_residues_ = setting; }

		void set_choose_random( bool const & setting ){ choose_random_ = setting; }
		bool choose_random() const{ return choose_random_; }

		void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }
		core::Size num_random_samples() const{ return num_random_samples_; }

		utility::vector1< core::id::TorsionID > const & which_torsions() const;

	private:

		StepWiseProteinCCD_CloserOP ccd_closer_;
		rotamer_sampler::RotamerSizedOP sampler_;

		utility::vector1< utility::vector1< core::Real > > main_chain_torsion_sets_;

		bool choose_random_;
		core::Size num_random_samples_;
		core::Size max_ntries_;

		core::Size ccd_close_res_;
		utility::vector1< core::Size > moving_residues_;

	};

} //loop_close
} //protein
} //sampling
} //stepwise
} //protocols

#endif
